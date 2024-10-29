date <- "240130MoreNoise"

source("../../src/1_generate_time.R")

## ----
## Set parameters
## ----
timing_path <- "timing/"
perf_score_path <- "perf_scores/"

# fs_bool = F for simulation with date later than 231027
generate_time_SB_silico(timing_path, perf_score_path, date)

## ----
## Verif values
## ----
vals = readRDS(paste0("perf_scores/",date,"_time.rds"))
range(vals$values, na.rm=T)
