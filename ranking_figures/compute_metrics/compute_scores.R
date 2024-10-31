source("generic_functions_metrics.R")

## ----
## Set parameters, change with your own
## ----
date <- "241025"
silico_path <- "../../data/simulations/"
vitro_path <- "../../data/invitro/"
vivo_path <- "../../data/invivo/"
deconv_path <- "../../deconvolution/results/prediction/"
timing_path <- "../../deconvolution/results/timing/"
score_path <- "scores/"
score_methods <- c("rmse", "mae", "pearson")
best2_deconv_methods <- data.frame(DeconvTool=c('rlr','nnls','ICA','NMF','rlr','WISP','debCAM','PREDE'),
                                   FS='none',
                                   Block=c(rep('dnam',4),rep('rna',4)),
                                   Supervised=rep(c(rep('sup',2),rep('unsup',2)),2))

## ----
## Scores
## ----
score_SB_silico(silico_path, deconv_path, score_path,
                  date, score_methods)
score_SB_consensus(silico_path, deconv_path, score_path,
                date, best2_deconv_methods, score_methods)
score_SB_invitro(vitro_path, deconv_path, score_path,
                   score_methods)
score_SB_invivo(vivo_path, deconv_path, score_path,
                   score_methods)

## ----
## Time
## ----
time_SB_silico(timing_path, score_path, date)
time_SB_consensus(silico_path, timing_path, score_path, date, best2_deconv_methods)
time_SB_invitro(timing_path, score_path)
time_SB_invivo(timing_path, score_path)
