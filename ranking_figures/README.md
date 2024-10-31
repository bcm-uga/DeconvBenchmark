## Contents
You will find here all the scripts to:
- compute the different metrics (RMSE, MAE, Pearson correlations, time) and store them in data frames: folder compute_metrics
- reproduce the different figures of the paper: folder figures

## Include new metrics
- modify generic_functions_metrics.R and compute_scores.R accordingly
- don't forget to set your own paths/parameters in compute_scores.R

## Modify the ranking process
The ranking process is computed for each figure with the script figures/generic_functions/ranking_process.R. To adapt the ranking process. one should modify:
- figures/generic_functions/ranking_process.R, where are all the functions needed to compute an overall score (normalization / aggregation procedure)
- figures/generic_functions/ranking_pval.R, where are the functions to do the permutation test