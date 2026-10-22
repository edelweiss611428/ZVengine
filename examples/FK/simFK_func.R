library("R6")
library("rstan")

compiled_model = stan_model(file = "./stan/FK.stan")
model_data     = rjson::fromJSON(file = "./data/FK.json")

initFK = function() {
  list(
    CL = exp(rnorm(1, log(10), 0.2)),
    Q = exp(rnorm(1, log(20), 0.2)),
    V1 = exp(rnorm(1, log(70), 0.2)),
    V2 = exp(rnorm(1, log(70), 0.2)),
    ka = exp(rnorm(1, log(1), 0.2)),
    sigma = runif(1, 0.5, 2),
    alpha = exp(rnorm(1, log(2e-3), 0.2)),
    mtt = exp(rnorm(1, log(125), 0.2)),
    circ0 = exp(rnorm(1, 5, 0.2)),
    gamma = exp(rnorm(1, 0.17, 0.2)),
    sigmaNeut = runif(1, 0.5, 2)
  )
}

simFK = function(i,
                 S = 100,
                 model = compiled_model,
                 data = model_data,
                 Qmax = 5,
                 K_set = c(25, 50),
                 n_warmup = 1000) {
  set.seed(i)
  
  t1 = proc.time()[3]
  fit = sampling(
    model,
    data = data,
    iter = n_warmup + S,
    warmup = n_warmup,
    chains = 1,
    thin = 1,
    control = list(adapt_delta = 0.95),
    init = initFK
  )
  t_mcmc = proc.time()[3] - t1
  
  ext = rstan::extract(fit)
  params = cbind(
    CL = ext$CL,
    Q = ext$Q,
    V1 = ext$V1,
    V2 = ext$V2,
    ka = ext$ka,
    sigma = ext$sigma,
    mtt = ext$mtt,
    circ0 = ext$circ0,
    alpha = ext$alpha,
    gamma = ext$gamma,
    sigmaNeut = ext$sigmaNeut
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
  save(results, file = sprintf("./FK/FK_S%d_%d.RData", S, i))
  invisible(results)
}
