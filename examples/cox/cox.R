library("rstan")
library(spatstat)

cox_stan_model = rstan::stan_model(file = "./stan/cox.stan")

make_cox_prior = function(ngrid) {
  dimension = ngrid^2
  
  parameter_sigmasq = 1.91
  parameter_mu = log(126) - 0.5 * parameter_sigmasq
  parameter_beta = 1 / 33
  parameter_area = 1 / dimension
  
  prior_mean = rep(parameter_mu, dimension)
  
  coords = cbind(rep(1:ngrid, each = ngrid), rep(1:ngrid, times = ngrid))
  
  dist_mat = as.matrix(dist(coords))
  
  prior_cov = parameter_sigmasq *
    exp(-dist_mat /
          (ngrid * parameter_beta))
  
  prior_cov_chol = t(chol(prior_cov))
  
  list(
    dimension = dimension,
    parameter_sigmasq = parameter_sigmasq,
    parameter_mu = parameter_mu,
    parameter_beta = parameter_beta,
    parameter_area = parameter_area,
    prior_mean = prior_mean,
    prior_cov = prior_cov,
    prior_cov_chol = prior_cov_chol
  )
}

make_cox_data = function(ngrid) {
  data(finpines)
  
  data_x = (finpines$x + 5) / 10
  data_y = (finpines$y + 8) / 10
  
  grid = seq(from = 0,
             to = 1,
             length.out = ngrid + 1)
  
  dimension = ngrid^2
  
  data_counts = rep(0, dimension)
  
  for (i in 1:ngrid) {
    for (j in 1:ngrid) {
      logical_y =
        (data_x > grid[i]) *
        (data_x < grid[i + 1])
      
      logical_x =
        (data_y > grid[j]) *
        (data_y < grid[j + 1])
      
      data_counts[(i - 1) * ngrid + j] = sum(logical_y * logical_x)
    }
  }
  
  list(
    data_x = data_x,
    data_y = data_y,
    grid = grid,
    data_counts = data_counts
  )
}

fit_cox = function(ngrid = 4,
                   chains = 1,
                   iter = 2000,
                   warmup = 1000,
                   seed = 1,
                   refresh = 100,
                   stan_model = cox_stan_model,
                   thin = 1) {
  obs = make_cox_data(ngrid = ngrid)
  
  prior = make_cox_prior(ngrid = ngrid)
  
  D = prior$dimension
  
  stan_data = list(
    D = D,
    data_counts = obs$data_counts,
    prior_mean = prior$prior_mean,
    L_prior_cov = prior$prior_cov_chol,
    area = prior$parameter_area
  )
  
  
  
  fitted_model = rstan::sampling(
    object = stan_model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    refresh = refresh
  )
  
  extracted_params = rstan::extract(fitted_model)
  
  params = extracted_params$x
  
  skel = get_inits(fitted_model)[[1]]
  
  grad = t(apply(
    params,
    MARGIN = 1,
    FUN = function(theta) {
      grad_log_prob(
        fitted_model,
        upars = unconstrain_pars(fitted_model, relist(theta, skel)),
        adjust_transform = TRUE
      )
    }
  ))
  
  unconstrained_pars = t(apply(
    params,
    MARGIN = 1,
    FUN = function(theta) {
      unconstrain_pars(fitted_model, relist(theta, skel))
    }
  ))
  
  list(
    data_x = obs$data_x,
    data_y = obs$data_y,
    data_counts = obs$data_counts,
    grid = obs$grid,
    prior_mean = prior$prior_mean,
    prior_cov = prior$prior_cov,
    prior_cov_chol = prior$prior_cov_chol,
    params = params,
    unconstrained_params = unconstrained_pars,
    grad = grad,
    ngrid = ngrid,
    dimension = D
  )
}
