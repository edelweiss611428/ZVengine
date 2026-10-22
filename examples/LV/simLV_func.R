library("R6")
library("rstan")

compiled_model = stan_model(file = "./stan/lv.stan")
model_data     = rjson::fromJSON(file = "./data/lv.json")

simLV = function(i,
                 S = 100,
                 model = compiled_model,
                 data = model_data,
                 Qmax = 5,
                 K_set = c(25, 50)) {
  set.seed(i)
  
  t1 = proc.time()[3]
  fit = sampling(
    model,
    data = data,
    iter = 1000 + S,
    warmup = 1000,
    chains = 1,
    thin = 1,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  t_mcmc = proc.time()[3] - t1
  
  ext = rstan::extract(fit)
  params = cbind(
    theta = ext$theta,
    z_init = ext$z_init,
    sigma = ext$sigma
  )
  skel = get_inits(fit)[[1]]
  upars = t(apply(params, 1, function(p)
    unconstrain_pars(fit, relist(p, skel))))
  grad = t(apply(upars, 1, function(u)
    grad_log_prob(fit, u, adjust_transform = TRUE)))
  
  d = ncol(upars)
  Q_ok = Filter(function(q)
    choose(q + d, d) < S, 1:Qmax)
  
  estimates = list(vanilla = colMeans(params))
  runtimes = list(vanilla = t_mcmc)
  
  run = function(name, expr) {
    t1 = proc.time()[3]
    val = as.vector(expr)
    runtimes[[name]] <<- proc.time()[3] - t1 + t_mcmc
    estimates[[name]] <<- val
  }
  skip = function(name) {
    runtimes[[name]] <<- NA
    estimates[[name]] <<- NA
  }
  
  for (q in 1:Qmax) {
    nm = paste0("ZV", q)
    if (q %in% Q_ok) {
      run(nm,
          ZVCV::zvcv(
            params,
            upars,
            grad,
            options = list(polyorder = q, regul_reg = FALSE)
          )$expectation)
    } else
      skip(nm)
    run(paste0("rZV", q),
        ZVCV::zvcv(
          params,
          upars,
          grad,
          options = list(
            polyorder = q,
            regul_reg = TRUE,
            alpha_elnet = 0
          )
        )$expectation)
    run(paste0("lZV", q),
        ZVCV::zvcv(
          params,
          upars,
          grad,
          options = list(
            polyorder = q,
            regul_reg = TRUE,
            alpha_elnet = 1
          )
        )$expectation)
  }
  
  run("CF", ZVCV::CF(params, upars, grad)$expectation)
  
  for (q in 1:2) {
    nm = paste0("SECF", q)
    if (q %in% Q_ok) {
      run(nm,
          ZVCV::SECF(
            params,
            upars,
            grad,
            polyorder = q,
            steinOrder = 2
          )$expectation)
    } else
      skip(nm)
  }
  
  Qbase = if (2 %in% Q_ok)
    1:2
  else
    1
  Qvary = (max(Qbase) + 1):Qmax
  
  engines = list(
    avg = SimpleAvgEngine$new(),
    qpp = QPCovarianceEngine$new(long_only = TRUE),
    qpn = QPCovarianceEngine$new(long_only = FALSE)
  )
  ens = function(opt, K, eng) {
    ZVengine$new(RandCand$new(
      Qbase = Qbase,
      Qvary = Qvary,
      options = opt,
      n = K
    ),
    eng)$fit(params, upars, grad)$ensemble_result$estimate
  }
  
  run("strat_1", ens("stratifiedEq", 1, engines$avg))
  
  opts = c(strat = "stratifiedEq", unif = "uniform")
  for (o in names(opts))
    for (K in K_set)
      for (e in names(engines)) {
        run(sprintf("%s_%s_%d", o, e, K), ens(opts[[o]], K, engines[[e]]))
      }
  
  estimates = as.data.frame(estimates, row.names = colnames(params))
  results = list(
    samples = list(
      unconstrain_pars = upars,
      params = params,
      grad = grad
    ),
    estimates = estimates,
    runtimes = runtimes
  )
  save(results, file = sprintf("./LV/LV_S%d_%d.RData", S, i))
  invisible(results)
}
