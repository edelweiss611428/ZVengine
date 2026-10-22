library(future)
library(future.apply)
library("parallel")

NRUN     = 100
NGRID    = 8
NCELL    = NGRID^2
S_VALUES = c(300, 1000)
NWORKERS = max(1, future::availableCores() - 1)
COX_DIR  = "cox"

method_names = c(
  "ZV1",
  "lZV1",
  "lZV2",
  "rZV1",
  "rZV2",
  "CF",
  "SECF1",
  "SECF2",
  "Anchor_k1",
  "Anchor_kmax",
  "Ours_SAvg",
  "Ours_QPplus"
)
all_names = c("vanilla", method_names)

source("examples/cox/cox.R")
source("ZVCV_engine.R")

library("spatstat")

rep_path = function(i, cox_dir = COX_DIR) {
  file.path(cox_dir, sprintf("rep_%04d.RData", i))
}

grid_nearest = function(index, n, k) {
  row = (index - 1) %/% n + 1
  col = (index - 1) %% n + 1
  grid_row = rep(1:n, each = n)
  grid_col = rep(1:n, times = n)
  distance = sqrt((grid_row - row)^2 + (grid_col - col)^2)
  grid_index = (grid_row - 1) * n + grid_col
  ord = order(distance, grid_index)
  grid_index[ord][1:min(k, n^2)]
}

max_feasible_k = function(S) {
  min(9, max(which(choose(0:20 + 4, 4) < S)) - 1)
}

# ------ Gold standard ------
get_gold_standard = function(ngrid = NGRID, cox_dir = COX_DIR) {
  gold_file = "examples/cox/gold_cox.csv"
  
  if (file.exists(gold_file)) {
    gold_standard = read.csv(gold_file)$gold
  } else {
    cat("Running gold-standard chain (S = 1e6) ...\n")
    gold_fit = fit_cox(
      warmup = 1e5,
      iter = 1e5 + 1e5,
      seed = 0,
      chains = 10,
      ngrid = ngrid
    )
    gold_standard = colMeans(gold_fit$params)
    write.csv(data.frame(parameter = paste0("x", seq_along(gold_standard)), gold = gold_standard),
              gold_file,
              row.names = FALSE)
  }
  gold_standard
}

# ------ Methods ------
run_cox_methods_for_S = function(S,
                                 params,
                                 usamples,
                                 grad,
                                 vanilla_time,
                                 NGRID) {
  NCELL = NGRID^2
  i = 1
  
  estimates_S = setNames(lapply(all_names, function(m)
    matrix(0, 1, NCELL)), all_names)
  runtimes_S  = setNames(lapply(all_names, function(m)
    numeric(1)), all_names)
  
  runtimes_S$vanilla[i]    = vanilla_time
  estimates_S$vanilla[i, ] = colMeans(params)
  
  t1 = proc.time()[3]
  estimates_S$ZV1[i, ] = as.vector(ZVCV::zvcv(
    params,
    usamples,
    grad,
    options = list(polyorder = 1, regul_reg = FALSE)
  )$expectation)
  t2 = proc.time()[3]
  runtimes_S$ZV1[i] = vanilla_time + (t2 - t1)
  
  for (q in 1:2) {
    t1 = proc.time()[3]
    estimates_S[[paste0("lZV", q)]][i, ] =
      as.vector(ZVCV::zvcv(
        params,
        usamples,
        grad,
        options = list(
          polyorder = q,
          regul_reg = TRUE,
          alpha_elnet = 1
        )
      )$expectation)
    t2 = proc.time()[3]
    runtimes_S[[paste0("lZV", q)]][i] = vanilla_time + (t2 - t1)
    
    t1 = proc.time()[3]
    estimates_S[[paste0("rZV", q)]][i, ] =
      as.vector(ZVCV::zvcv(
        params,
        usamples,
        grad,
        options = list(
          polyorder = q,
          regul_reg = TRUE,
          alpha_elnet = 0
        )
      )$expectation)
    t2 = proc.time()[3]
    runtimes_S[[paste0("rZV", q)]][i] = vanilla_time + (t2 - t1)
  }
  
  t1 = proc.time()[3]
  estimates_S$CF[i, ] = as.vector(ZVCV::CF(params, usamples, grad)$expectation)
  t2 = proc.time()[3]
  runtimes_S$CF[i] = vanilla_time + (t2 - t1)
  
  t1 = proc.time()[3]
  estimates_S$SECF1[i, ] = as.vector(ZVCV::SECF(
    params,
    usamples,
    grad,
    polyorder = 1,
    steinOrder = 2
  )$expectation)
  t2 = proc.time()[3]
  runtimes_S$SECF1[i] = vanilla_time + (t2 - t1)
  
  t1 = proc.time()[3]
  estimates_S$SECF2[i, ] = as.vector(ZVCV::SECF(
    params,
    usamples,
    grad,
    polyorder = 2,
    steinOrder = 2
  )$expectation)
  t2 = proc.time()[3]
  runtimes_S$SECF2[i] = vanilla_time + (t2 - t1)
  
  ZZ = max_feasible_k(S)
  
  t1 = proc.time()[3]
  for (j in 1:NCELL) {
    anchor_spec = aprioriZV$new(
      degree = 1:4,
      fixed_dims = j,
      n_random_dims = 0,
      max_dims = 1,
      n = 1
    )
    eng = ZVengine$new(specs = list(anchor_spec),
                       combination_engine = SimpleAvgEngine$new())
    eng$fit(y = params[, j],
            samples = usamples,
            grad = grad)
    estimates_S$Anchor_k1[i, j] = eng$ensemble_result$estimate
  }
  t2 = proc.time()[3]
  runtimes_S$Anchor_k1[i] = vanilla_time + (t2 - t1)
  
  t1 = proc.time()[3]
  for (j in 1:NCELL) {
    nk = grid_nearest(index = j, n = NGRID, k = ZZ)
    anchor_spec = aprioriZV$new(
      degree = 1:4,
      fixed_dims = nk,
      n_random_dims = 0,
      max_dims = length(nk),
      n = 1
    )
    eng = ZVengine$new(specs = list(anchor_spec),
                       combination_engine = SimpleAvgEngine$new())
    eng$fit(y = params[, j],
            samples = usamples,
            grad = grad)
    estimates_S$Anchor_kmax[i, j] = eng$ensemble_result$estimate
  }
  t2 = proc.time()[3]
  runtimes_S$Anchor_kmax[i] = vanilla_time + (t2 - t1)
  
  fit_time_ours = 0
  combine_time_QPplus = 0
  
  for (j in 1:NCELL) {
    specs = lapply(1:ZZ, function(k) {
      nk = grid_nearest(index = j, n = NGRID, k = k)
      aprioriZV$new(
        degree = 1:4,
        fixed_dims = nk,
        n_random_dims = 0,
        max_dims = length(nk),
        n = 1
      )
    })
    
    t1 = proc.time()[3]
    ours = ZVengine$new(specs = specs, combination_engine = SimpleAvgEngine$new())
    ours$fit(y = params[, j],
             samples = usamples,
             grad = grad)
    t2 = proc.time()[3]
    fit_time_ours = fit_time_ours + (t2 - t1)
    estimates_S$Ours_SAvg[i, j] = ours$ensemble_result$estimate
    
    t1 = proc.time()[3]
    est_QPplus = QPCovarianceEngine$new(long_only = TRUE)$combine(ours$estimates, ours$loo_residuals)
    t2 = proc.time()[3]
    combine_time_QPplus = combine_time_QPplus + (t2 - t1)
    estimates_S$Ours_QPplus[i, j] = est_QPplus$estimate
  }
  
  runtimes_S$Ours_SAvg[i]   = vanilla_time + fit_time_ours
  runtimes_S$Ours_QPplus[i] = vanilla_time + fit_time_ours + combine_time_QPplus
  
  list(estimates = estimates_S, runtimes = runtimes_S)
}

# ------ One replicate ------
run_one_rep = function(i, S_values, NGRID_, cox_dir) {
  f = file.path(cox_dir, sprintf("rep_%04d.RData", i))
  if (file.exists(f)) {
    ok = tryCatch({
      load(f)
      TRUE
    }, error = function(e)
      FALSE)
    if (ok && exists("rep_out"))
      return(rep_out)
  }
  
  source("examples/cox/cox.R")
  source("ZVCV_engine.R")
  
  S_max = max(S_values)
  
  t1 = proc.time()[3]
  fit = fit_cox(
    warmup = 1000,
    iter = 1000 + S_max,
    seed = i,
    chains = 1,
    ngrid = NGRID_
  )
  t2 = proc.time()[3]
  full_time = t2 - t1
  
  grad_full     = fit$grad
  usamples_full = fit$unconstrained_params
  params_full   = fit$params
  
  rep_out = setNames(lapply(S_values, function(S) {
    idx = seq_len(S)
    run_cox_methods_for_S(
      S = S,
      params   = params_full[idx, , drop = FALSE],
      usamples = usamples_full[idx, , drop = FALSE],
      grad     = grad_full[idx, , drop = FALSE],
      vanilla_time = full_time * (S / S_max),
      NGRID = NGRID_
    )
  }), paste0("S", S_values))
  
  tmp = paste0(f, ".tmp")
  save(rep_out, file = tmp)
  file.rename(tmp, f)
  
  rep_out
}

# ------ Parallel replicates ------
run_cox_multi_S = function(S_values = S_VALUES,
                           NRUN_ = NRUN,
                           NGRID_ = NGRID,
                           workers = NWORKERS,
                           cox_dir = COX_DIR) {
  dir.create(cox_dir, showWarnings = FALSE, recursive = TRUE)
  
  done = sum(file.exists(vapply(
    1:NRUN_, rep_path, character(1), cox_dir = cox_dir
  )))
  if (done > 0)
    cat("Found",
        done,
        "completed replicates in",
        cox_dir,
        "- reusing them.\n")
  
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  
  cat("Running",
      NRUN_ - done,
      "replicates on",
      workers,
      "workers ...\n")
  t0 = proc.time()[3]
  
  reps = future_lapply(
    1:NRUN_,
    run_one_rep,
    S_values = S_values,
    NGRID_ = NGRID_,
    cox_dir = cox_dir,
    future.seed = TRUE
  )
  
  cat(sprintf("Done in %.1f min\n", (proc.time()[3] - t0) / 60))
  
  estimates = setNames(lapply(paste0("S", S_values), function(Sname)
    setNames(lapply(all_names, function(m)
      do.call(rbind, lapply(reps, function(r)
        r[[Sname]]$estimates[[m]]))), all_names)), paste0("S", S_values))
  
  runtimes = setNames(lapply(paste0("S", S_values), function(Sname)
    setNames(lapply(all_names, function(m)
      vapply(reps, function(r)
        r[[Sname]]$runtimes[[m]], numeric(1))), all_names)), paste0("S", S_values))
  
  save(estimates,
       runtimes,
       file = file.path(cox_dir, "cox_multi_S_progress.RData"))
  
  list(
    estimates = estimates,
    runtimes = runtimes,
    method_names = method_names,
    S_values = S_values
  )
}

# ------ Summarise ------
summarise_cox = function(run_out, gold_standard) {
  S_values     = run_out$S_values
  method_names = run_out$method_names
  
  long = do.call(rbind, lapply(S_values, function(S) {
    Sname = paste0("S", S)
    est = run_out$estimates[[Sname]]
    rt  = run_out$runtimes[[Sname]]
    
    mse_vanilla = colMeans(sweep(est$vanilla, 2, gold_standard, "-")^2)
    
    rows = lapply(method_names, function(m) {
      mse_m = colMeans(sweep(est[[m]], 2, gold_standard, "-")^2)
      SE_m  = mean(mse_vanilla / mse_m)
      CE_m  = mean(rt$vanilla / rt[[m]])
      data.frame(
        method = m,
        S = S,
        SE = SE_m,
        CE = CE_m,
        Overall = SE_m * CE_m
      )
    })
    do.call(rbind, rows)
  }))
  rownames(long) = NULL
  
  wide = reshape(long[, c("method", "S", "Overall")],
                 idvar = "method",
                 timevar = "S",
                 direction = "wide")
  names(wide) = sub("^Overall\\.", "S", names(wide))
  wide = wide[match(method_names, wide$method), ]
  
  list(long = long, wide = wide)
}

# ------ LaTeX table ------
tex_label = c(
  ZV1 = "$\\text{ZV}_1$",
  lZV1 = "$\\ell\\text{-}\\text{ZV}_1$",
  lZV2 = "$\\ell\\text{-}\\text{ZV}_2$",
  rZV1 = "$r\\text{-}\\text{ZV}_1$",
  rZV2 = "$r\\text{-}\\text{ZV}_2$",
  CF = "$\\text{CF}$",
  SECF1 = "$\\text{SECF}_1$",
  SECF2 = "$\\text{SECF}_2$",
  Anchor_k1 = "$\\text{Anchor}_{k=1}$",
  Anchor_kmax = "$\\text{Anchor}_{k=\\max}$",
  Ours_SAvg = "$\\mathrm{Ours}_{\\mathrm{SAvg}}$",
  Ours_QPplus = "$\\mathrm{Ours}_{\\mathrm{QP}_{+}}$"
)

make_latex_table_cox = function(summary,
                                method_names,
                                metric = "Overall",
                                digits = 2,
                                S_values = NULL) {
  long = summary$long
  if (is.null(S_values))
    S_values = sort(unique(long$S))
  
  lines = sapply(method_names, function(m) {
    vals = sapply(S_values, function(S) {
      row = long[long$method == m & long$S == S, ]
      if (nrow(row) == 0)
        NA
      else
        row[[metric]]
    })
    vals_str = ifelse(is.na(vals),
                      "--",
                      formatC(vals, digits = digits, format = "f"))
    paste0(tex_label[m], " & ", paste(vals_str, collapse = " & "), " \\\\")
  })
  
  cat(paste(lines, collapse = "\n"), "\n\\bottomrule\n")
  invisible(lines)
}

# ------ Run ------
run_and_summarise_cox = function(S_values = S_VALUES,
                                 NRUN_ = NRUN,
                                 NGRID_ = NGRID,
                                 metric = "Overall",
                                 digits = 2,
                                 workers = NWORKERS,
                                 cox_dir = COX_DIR) {
  dir.create(cox_dir, showWarnings = FALSE, recursive = TRUE)
  gold_standard = get_gold_standard(ngrid = NGRID_, cox_dir = cox_dir)
  
  run_out = run_cox_multi_S(
    S_values = S_values,
    NRUN_ = NRUN_,
    NGRID_ = NGRID_,
    workers = workers,
    cox_dir = cox_dir
  )
  summary = summarise_cox(run_out, gold_standard)
  
  cat("\n=== Long table (method x S, SE / CE / Overall) ===\n")
  print(summary$long)
  
  cat("\n=== Wide table (method x S, Overall = SE * CE) ===\n")
  print(summary$wide)
  
  cat("\n=== LaTeX rows (metric =", metric, ") ===\n")
  make_latex_table_cox(
    summary,
    method_names = run_out$method_names,
    metric = metric,
    digits = digits,
    S_values = S_values
  )
  
  write.csv(summary$long,
            file.path(cox_dir, "cox_multi_S_summary_long.csv"),
            row.names = FALSE)
  write.csv(summary$wide,
            file.path(cox_dir, "cox_multi_S_summary_wide.csv"),
            row.names = FALSE)
  
  invisible(list(run_out = run_out, summary = summary))
}

res = run_and_summarise_cox(S_values = c(300, 1000), NRUN_ = 100)
