# Ensemble ZVCV engine

Run everything from the project root (open `nEE.Rproj`; create one if this does not exist).

```
ZVCV_engine.R, getXfast.cpp   ensemble ZVCV engine
stan/                         Stan models (lv, FK, cox)
data/                         Stan data (lv.json, FK.json)
examples/<LV|FK|cox>/         scripts and gold standards
```

## Lotka–Volterra / Friberg–Karlsson

Results: set `S` (100, 300, 1000, 3000) in the parallel script, then

```r
source("examples/LV/simLV100.R")     # reps 1:101, rep 32 fails -> 100 reps
source("examples/FK/simFK_100.R")    # reps 1:100
```

Per-replicate files go to `LV/` and `FK/`. Tables and efficiency CSVs:

```r
source("examples/LV/summarise_LV.R")
source("examples/FK/summarise_FK.R")
```

Gold standards: `examples/LV/gold_LV.csv`, `examples/FK/gold_FK.csv` (precomputed).

## Cox process

```r
source("examples/cox/cox_experiment.R")   # 100 reps, S = 300, 1000; summary + LaTeX rows
```

Gold standard: `examples/cox/gold_cox.csv` (precomputed); if absent, `get_gold_standard()` recomputes it from a long chain.
