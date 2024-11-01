## Contents
You will find here all the scripts to:
- compute the different metrics (RMSE, MAE, Pearson correlations, time) and store them in data frames: folder [compute_metrics](compute_metrics/),
- reproduce the different figures of the paper: folder [figures](figures/).

## Include new metrics
- modify [generic_functions_metrics.R](compute_metrics/generic_functions_metrics.R) and [compute_scores.R](compute_metrics/compute_scores.R) accordingly
- don't forget to set your own paths/parameters in [compute_scores.R](compute_metrics/compute_scores.R).

## Modify the ranking process
The ranking process is computed for each figure with the script [ranking_process.R](figures/generic_functions/ranking_process.R). To adapt the ranking process. one should modify:
- [ranking_process.R](figures/generic_functions/ranking_process.R), where all the functions needed to compute an overall score (normalization / aggregation procedure) are,
- [ranking_pval.R](figures/generic_functions/ranking_pval.R), where the functions to do the permutation test are.
