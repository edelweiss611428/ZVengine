gold = with(read.csv("examples/FK/gold_FK.csv"),
            setNames(gold, parameter))

dir    = "./FK"
S_set  = c(100, 300, 1000, 3000)
reps   = 1:100
outTex = "FK_tables.tex"

readRep = function(S, i) {
  f = sprintf("%s/FK_S%d_%d.RData", dir, S, i)
  if (!file.exists(f))
    return(NULL)
  e = new.env()
  load(f, envir = e)
  e$results
}

effS = function(S) {
  out = Filter(Negate(is.null), lapply(reps, readRep, S = S))
  if (!length(out))
    return(NULL)
  methods = colnames(out[[1]]$estimates)
  pars    = rownames(out[[1]]$estimates)
  g = if (all(names(gold) %in% pars))
    gold[pars]
  else {
    stopifnot(length(pars) == length(gold))
    gold
  }
  
  E = simplify2array(lapply(out, function(r)
    as.matrix(r$estimates)))
  R = t(sapply(out, function(r)
    vapply(r$runtimes[methods], function(x)
      as.numeric(x)[1], numeric(1))))
  colnames(R) = methods
  
  mse  = apply((E - array(g, dim(E)))^2, 1:2, mean)
  t_mc = mean(R[, "vanilla"])
  cost = pmax(colMeans(R) - t_mc, 0)
  se   = colMeans(mse[, "vanilla"] / mse)
  
  data.frame(
    method = methods,
    SE = se,
    OE = se * t_mc / (cost + t_mc),
    cost = cost,
    t_mc = t_mc,
    nrep = length(out),
    row.names = NULL
  )
}

tabs = Filter(Negate(is.null), setNames(lapply(S_set, effS), S_set))
S_ok = as.numeric(names(tabs))

fmt = function(x)
  ifelse(is.na(x), "--", formatC(x, format = "e", digits = 2))

emit = function(keys,
                labels,
                groups,
                ours,
                caption,
                label,
                width = 0.8) {
  se = sapply(tabs, function(t)
    t$SE[match(keys, t$method)])
  oe = sapply(tabs, function(t)
    t$OE[match(keys, t$method)])
  M  = cbind(matrix(se, length(keys)), matrix(oe, length(keys)))
  nS = length(S_ok)
  
  body = character()
  for (r in seq_along(keys)) {
    cells = fmt(M[r, ])
    for (j in seq_len(ncol(M)))
      if (!is.na(M[r, j]) && M[r, j] == max(M[, j], na.rm = TRUE))
        cells[j] = sprintf("\\textcolor{red}{\\textbf{%s}}", cells[j])
    lab = if (keys[r] %in% ours)
      sprintf("$\\bm{%s}$", labels[r])
    else
      sprintf("$%s$", labels[r])
    body = c(body, sprintf("%s & %s \\\\", lab, paste(cells, collapse = " & ")))
    if (r %in% groups)
      body = c(body, "\\midrule")
  }
  
  c(
    "\\begin{table}[h]",
    "    \\renewcommand{\\arraystretch}{1}",
    "    \\centering",
    sprintf("    \\caption{%s}", caption),
    sprintf("    \\label{%s}", label),
    sprintf("    \\resizebox{%s\\textwidth}{!}{", width),
    sprintf("\\begin{tabular}{c%s}", strrep("c", 2 * nS)),
    "\\toprule",
    sprintf(
      "& \\multicolumn{%d}{c}{Statistical efficiency} & \\multicolumn{%d}{c}{Overall efficiency} \\\\",
      nS,
      nS
    ),
    sprintf(
      "\\cmidrule(r){2-%d}\\cmidrule(l){%d-%d}",
      nS + 1,
      nS + 2,
      2 * nS + 1
    ),
    sprintf("Method & %s \\\\", paste(rep(
      sprintf("$S = %d$", S_ok), 2
    ), collapse = " & ")),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\end{table}",
    ""
  )
}

ens   = c(
  "strat_1",
  "strat_avg_25",
  "strat_qpp_25",
  "strat_qpn_25",
  "strat_avg_50",
  "strat_qpp_50",
  "strat_qpn_50"
)
keys1 = c(
  sprintf("ZV%d", 1:5),
  sprintf("lZV%d", 1:5),
  sprintf("rZV%d", 1:5),
  "CF",
  "SECF1",
  "SECF2",
  ens
)
lab1  = c(
  sprintf("\\text{ZV}_%d", 1:5),
  sprintf("\\ell\\text{-}\\text{ZV}_%d", 1:5),
  sprintf("r\\text{-}\\text{ZV}_%d", 1:5),
  "\\text{CF}",
  "\\text{SECF}_1",
  "\\text{SECF}_2",
  "\\mathrm{eZV}_{1}",
  "\\mathrm{eZV}_{25,\\mathrm{SAvg}}",
  "\\mathrm{eZV}_{25,\\mathrm{QP}_{+}}",
  "\\mathrm{eZV}_{25,\\mathrm{QP}_{\\pm}}",
  "\\mathrm{eZV}_{50,\\mathrm{SAvg}}",
  "\\mathrm{eZV}_{50,\\mathrm{QP}_{+}}",
  "\\mathrm{eZV}_{50,\\mathrm{QP}_{\\pm}}"
)
grp1  = cumsum(c(5, 5, 5, 3))

eng   = c(avg = "\\mathrm{SAvg}", qpp = "\\mathrm{QP}_{+}", qpn = "\\mathrm{QP}_{\\pm}")
comb  = expand.grid(
  e = names(eng),
  K = c(25, 50),
  o = c("strat", "unif"),
  stringsAsFactors = FALSE
)
keys2 = sprintf("%s_%s_%d", comb$o, comb$e, comb$K)
lab2  = sprintf("\\mathrm{eZV}^{\\mathrm{%s}}_{%d,%s}", comb$o, comb$K, eng[comb$e])
grp2  = nrow(comb) / 2

tex = c(
  emit(
    keys1,
    lab1,
    grp1,
    ens,
    "Results of Experiment 2: Friberg--Karlsson Model ($d = 11$)",
    "table:fk"
  ),
  emit(
    keys2,
    lab2,
    grp2,
    keys2,
    "Experiment 2: stratified versus uniform basis selection ($d = 11$)",
    "table:fk-strat-unif",
    width = 0.75
  )
)

writeLines(tex, outTex)
cat(tex, sep = "\n")

write.csv(do.call(rbind, lapply(names(tabs), function(s)
  transform(tabs[[s]], S = as.numeric(s)))), "FK_efficiency.csv", row.names = FALSE)
