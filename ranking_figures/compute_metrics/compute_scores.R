source("generic_functions_metrics.R")

## ----
## Set parameters
## ----
date <- "241025"
data_path <- "../../data/simulations/"
deconv_path <- "../../deconvolution/results/prediction/"
score_path <- "scores/"
score_methods <- c("rmse", "mae", "pearson")

## ----
## Scores
## ----
setting_SB_silico(data_path, deconv_path, score_path,
                  date, score_methods)
