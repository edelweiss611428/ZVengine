library(R6)
Rcpp::sourceCpp("./getXfast.cpp")

# ------ Helpers ------
getCp = function(d, Q) {
  if (!length(Q))
    return(matrix(0, d, 0))
  do.call(cbind, lapply(Q, function(q)
    as.matrix(partitions::compositions(q, d))))
}

# ------ Weak learners ------
WeakLearner = R6Class("WeakLearner", public = list(
  generate_features = function(samples, grad)
    stop("Not implemented.")
))

ZVCV = R6Class(
  "ZVCV",
  inherit = WeakLearner,
  public = list(
    degree = NULL,
    n_instances = NULL,
    
    initialize = function(degree = 1, n = 1) {
      self$degree = degree
      self$n_instances = n
    },
    
    generate_features = function(samples, grad) {
      Z = efficient_getX(samples, grad, getCp(ncol(samples), self$degree))
      replicate(self$n_instances, Z, simplify = FALSE)
    }
  )
)

RandCand = R6Class(
  "RandCand",
  inherit = WeakLearner,
  public = list(
    Qbase = NULL,
    Qvary = NULL,
    options = NULL,
    n_candidates = NULL,
    max_feats_override = NULL,
    n_instances = NULL,
    
    initialize = function(Qbase = 1:2,
                          Qvary = 3:5,
                          options = "stratifiedEq",
                          n = 20,
                          max_feats = NULL,
                          n_instances = 1) {
      self$Qbase = Qbase
      self$Qvary = Qvary
      self$options = match.arg(options, c("uniform", "stratifiedEq"))
      self$n_candidates = n
      self$max_feats_override = max_feats
      self$n_instances = as.integer(n_instances)
    },
    
    generate_features = function(samples, grad) {
      S = nrow(samples)
      d = ncol(samples)
      
      Z_base = efficient_getX(samples, grad, getCp(d, self$Qbase))
      M_max = if (is.null(self$max_feats_override))
        floor(min(0.8 * S, 25 * sqrt(S)))
      else
        self$max_feats_override
      M_rem = max(0, M_max - ncol(Z_base))
      
      vary_indices = list()
      Cp_list = list()
      offset = 0
      for (q in self$Qvary) {
        Cp_q = getCp(d, q)
        n_q = ncol(Cp_q)
        if (n_q > 0) {
          Cp_list[[length(Cp_list) + 1]] = Cp_q
          vary_indices[[as.character(q)]] = offset + seq_len(n_q)
          offset = offset + n_q
        }
      }
      total_vary = offset
      Z_vary = if (total_vary > 0)
        efficient_getX(samples, grad, do.call(cbind, Cp_list))
      else
        NULL
      
      quota = setNames(numeric(length(vary_indices)), names(vary_indices))
      if (self$options == "stratifiedEq" &&
          M_rem > 0 && total_vary > 0) {
        budget = min(M_rem, total_vary)
        active = names(vary_indices)
        counts = lengths(vary_indices)
        
        while (budget > 0 && length(active) > 0) {
          target = floor(budget / length(active))
          if (target == 0) {
            chosen = if (length(active) == 1)
              active
            else
              sample(active, budget)
            quota[chosen] = quota[chosen] + 1
            break
          }
          add = pmin(counts[active] - quota[active], target)
          quota[active] = quota[active] + add
          budget = budget - sum(add)
          active = active[counts[active] > quota[active]]
        }
      }
      
      feature_list = vector("list", self$n_candidates)
      for (k in seq_len(self$n_candidates)) {
        cols = integer(0)
        if (M_rem > 0 && total_vary > 0) {
          if (self$options == "uniform") {
            cols = sample.int(total_vary, min(M_rem, total_vary))
          } else {
            cols = unlist(lapply(names(quota)[quota > 0], function(deg) {
              pool = vary_indices[[deg]]
              if (quota[deg] == length(pool))
                pool
              else
                pool[sample.int(length(pool), quota[deg])]
            }), use.names = FALSE)
          }
        }
        feature_list[[k]] = if (length(cols))
          cbind(Z_base, Z_vary[, cols, drop = FALSE])
        else
          Z_base
      }
      
      if (self$n_instances > 1)
        feature_list = rep(feature_list, each = self$n_instances)
      feature_list
    }
  )
)

aprioriZV = R6Class(
  "aprioriZV",
  inherit = WeakLearner,
  public = list(
    degree = NULL,
    fixed_dims = NULL,
    n_random_dims = NULL,
    max_dims = NULL,
    n_candidates = NULL,
    
    initialize = function(degree = 1,
                          fixed_dims = NULL,
                          n_random_dims = 0,
                          max_dims = NULL,
                          n = 1) {
      self$degree = degree
      self$fixed_dims = if (is.null(fixed_dims))
        integer(0)
      else
        as.integer(fixed_dims)
      self$n_random_dims = as.integer(n_random_dims)
      self$max_dims = if (is.null(max_dims))
        NULL
      else
        as.integer(max_dims)
      self$n_candidates = as.integer(n)
    },
    
    generate_features = function(samples, grad) {
      d = ncol(samples)
      fixed = self$fixed_dims
      if (any(fixed < 1 |
              fixed > d))
        stop("All 'fixed_dims' must lie in 1..d.")
      
      avail = setdiff(seq_len(d), fixed)
      n_draw = min(self$n_random_dims, length(avail))
      
      feature_list = vector("list", self$n_candidates)
      for (k in seq_len(self$n_candidates)) {
        rand = if (n_draw == 0)
          integer(0)
        else if (length(avail) == 1)
          avail
        else
          sample(avail, n_draw)
        
        dims = sort(unique(c(fixed, rand)))
        if (!is.null(self$max_dims) &&
            length(dims) > self$max_dims) {
          dims = if (length(fixed) >= self$max_dims)
            fixed[seq_len(self$max_dims)]
          else
            c(fixed, rand[seq_len(self$max_dims - length(fixed))])
          dims = sort(dims)
        }
        
        feature_list[[k]] = if (!length(dims)) {
          matrix(0, nrow(samples), 0)
        } else {
          efficient_getX(samples[, dims, drop = FALSE], grad[, dims, drop = FALSE], getCp(length(dims), self$degree))
        }
      }
      feature_list
    }
  )
)

# ------ Combination engines ------
CombinationEngine = R6Class("CombinationEngine", public = list(
  combine = function(estimates, loo_residuals)
    stop("Not implemented.")
))

combine_per_target = function(estimates, loo_residuals, weight_fn) {
  estimates = as.matrix(estimates)
  K = nrow(estimates)
  M = ncol(estimates)
  if (length(dim(loo_residuals)) == 2)
    dim(loo_residuals) = c(dim(loo_residuals), 1)
  
  W = matrix(0, K, M)
  for (m in seq_len(M)) {
    R = matrix(loo_residuals[, , m], ncol = K)
    v = apply(R, 2, var)
    W[, m] = if (any(v < 1e-12))
      replace(numeric(K), which.min(v), 1)
    else
      weight_fn(cov(R))
  }
  list(estimate = colSums(W * estimates), weights = W)
}

gmv_weights = function(Sigma) {
  K = nrow(Sigma)
  x = solve((Sigma + t(Sigma)) / 2 + diag(1e-8, K), rep(1, K))
  as.vector(x / sum(x))
}

qp_weights = function(Sigma) {
  K = nrow(Sigma)
  quadprog::solve.QP(
    Dmat = Sigma + diag(1e-8, K),
    dvec = rep(0, K),
    Amat = cbind(1, diag(K)),
    bvec = c(1, rep(0, K)),
    meq = 1
  )$solution
}

SimpleAvgEngine = R6Class("SimpleAvgEngine",
                          inherit = CombinationEngine,
                          public = list(
                            combine = function(estimates, loo_residuals) {
                              estimates = as.matrix(estimates)
                              list(estimate = colMeans(estimates),
                                   weights = matrix(1 / nrow(estimates), nrow(estimates), ncol(estimates)))
                            }
                          ))

ClassicalGMVPEngine = R6Class(
  "ClassicalGMVPEngine",
  inherit = CombinationEngine,
  public = list(
    combine = function(estimates, loo_residuals)
      combine_per_target(estimates, loo_residuals, gmv_weights)
  )
)

QPCovarianceEngine = R6Class(
  "QPCovarianceEngine",
  inherit = CombinationEngine,
  public = list(
    long_only = TRUE,
    
    initialize = function(long_only = TRUE)
      self$long_only = long_only,
    
    combine = function(estimates, loo_residuals)
      combine_per_target(estimates, loo_residuals, if (self$long_only)
        qp_weights
        else
          gmv_weights)
  )
)

# ------ ZVengine ------
ZVengine = R6Class(
  "ZVengine",
  public = list(
    specs = NULL,
    combination_engine = NULL,
    estimates = NULL,
    loo_residuals = NULL,
    ensemble_result = NULL,
    
    initialize = function(specs = list(),
                          combination_engine = SimpleAvgEngine$new()) {
      if (!is.list(specs))
        specs = list(specs)
      self$specs = specs
      self$combination_engine = combination_engine
    },
    
    fit = function(y, samples, grad) {
      Y = as.matrix(y)
      S = nrow(Y)
      M = ncol(Y)
      
      designs = unlist(lapply(self$specs, function(s)
        s$generate_features(samples, grad)),
        recursive = FALSE)
      K = length(designs)
      
      self$estimates = matrix(0, K, M)
      self$loo_residuals = array(0, c(S, K, M))
      for (k in seq_len(K)) {
        fit_k = fast_ols_loo_cpp(Y, designs[[k]])
        self$estimates[k, ] = as.vector(fit_k$intercept)
        self$loo_residuals[, k, ] = fit_k$loo_residuals
      }
      
      self$ensemble_result = self$combination_engine$combine(self$estimates, self$loo_residuals)
      invisible(self)
    }
  )
)
