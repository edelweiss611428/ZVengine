library("doParallel")
library("foreach")
library("future")
library("rstan")
library("R6")
library("ZVCV")
S = 100
reps = 1:101

compiled_model = stan_model(file = "./stan/lv.stan")
model_data = rjson::fromJSON(file = "./data/lv.json")

n_core = 15
dir.create("LV", showWarnings = FALSE)

cl = makeCluster(n_core)
registerDoParallel(cl)
plan(multisession)

t0 = proc.time()[3]

out = foreach(i = reps) %dopar% {
  tryCatch({
    library(R6)
    library(rstan)
    library(ZVCV)
    
    rstan_options(auto_write = TRUE)
    options(mc.cores = 1)
    source("ZVCV_engine.R")
    source("examples/LV/simLV_func.R")
    
    simLV(i, S = S)
    
    "ok"
  }, error = function(e) {
    conditionMessage(e)
  })
}

stopCluster(cl)

status = setNames(unlist(out), reps)

cat(sprintf(
  "Done in %.1f min on %d workers; %d/%d ok\n",
  (proc.time()[3] - t0) / 60,
  n_core,
  sum(status == "ok"),
  length(reps)
))

print(status[status != "ok"])
