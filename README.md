# A robust workflow to benchmark deconvolution of multi-omic data

This repository contains all the files needed to perform deconvolution on example datasets, along with the ranking and the figures of the paper.
Our pipeline can be easily extended to include and evaluate novel methods, as well as other datasets.

## How to make in silico data: folder [data](data/)

Shortly, we made in silico data using reference profiles of pure cell types from different tissues convoluted with proportions generated based on a Dirichlet distribution. The scripts are in [simulation_scripts](data/simulation_scripts/), and the reference profiles can be downloaded from Zenodo (DOI 10.5281/zenodo.14024479).
The folder [data](data/) has the following architecture:
```
.
├── simulation_scripts
│   ├── generate_simu_DATA1.R
│   └── generate_simu_DATA2.R
├── references
│   ├── DATA1.rds
│   ├── DATA2_dnam.rds
│   └── DATA2_rna.rds
├── simulations
│   ├── dnam
│   │   ├── 240101_DATA1_sim01.rds
│   │   ├── 240101_DATA1_sim02.rds
│   │   ├── 240101_DATA2_sim01.rds
│   │   └── 240101_DATA2_sim02.rds
│   ├── rna
│   │   ├── 240101_DATA2_sim01.rds
│   │   └── 240101_DATA2_sim02.rds
...
```

This is also where we stored the *in vitro* and *in vivo* datasets, with the syntax 
```data/invitro/DATA_A.rds``` for the reference matrix and ```data/invitro/DATA_D_OMIC.rds``` for the methylation/expression matrix/. Please refer to the table [here](data/README.md) for the instructions on where to download *in vitro* and *in vivo* data.

## How to run the deconvolution methods: folder [deconvolution](deconvolution/)

This pipeline uses an Apptainer container and Snakemake. Instructions on how to use these tools can be found [here](deconvolution/README.md) and [here](deconvolution/container.md).

Briefly, this folder contains the scripts to perform the deconvolution pipeline. There is one script per setting (class of the method/omic type). Those scripts can be modified to include new methods and/or new datasets (cf [README](deconvolution/README.md)).

Snakemake will run all methods for all omics. The [Snakefile](deconvolution/Snakefile) is self-explanatory and can be modified to include new methods/datasets. In general, you can refer to the [README](deconvolution/README.md) to know how to test new methods/datasets.

Results of the deconvolution, *i.e.* estimation of the proportion matrix along with elapsed time will be stored in ```deconvolution/results/prediction/OMIC/CLASS/``` for the proportion matrix and ```deconvolution/results/timing/OMIC/CLASS/``` for the time elapsed with the syntax ```240101_DATA1_Apred_FS_METHOD_sim01.rds``` / ```240101_DATA1_timing_FS_METHOD_sim01.rds``` (FS being the feature selection strategy and METHOD the deconvolution algorithm).

## How to do the ranking and reproduce the figures of the paper: folder [ranking_figures](ranking_figures/)

(a) First, you can compute the different metrics (in our case, RMSE, MAE and Pearson correlation coefficients): just run the script [compute_scores.R](compute_metrics/compute_scores.R) and the scores will be stored in ```compute_metrics/scores/```: one file for the time (```..._time.rds```) and one file for the other metrics (```..._scores.rds```)

(b) The different figures of the paper can then be reproduced.

## Session Info
```R
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 14.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

loaded via a namespace (and not attached):
 [1] rstudioapi_0.16.0     knitr_1.46           
 [3] cluster_2.1.4         BiocGenerics_0.44.0  
 [5] IRanges_2.32.0        doParallel_1.0.17    
 [7] colorspace_2.1-0      clue_0.3-65          
 [9] rjson_0.2.21          rlang_1.1.3          
[11] fastmap_1.1.1         foreach_1.5.2        
[13] tools_4.2.0           parallel_4.2.0       
[15] grid_4.2.0            circlize_0.4.15      
[17] xfun_0.43             png_0.1-8            
[19] cli_3.6.2             htmltools_0.5.8.1    
[21] iterators_1.0.14      matrixStats_1.2.0    
[23] yaml_2.3.8            digest_0.6.35        
[25] crayon_1.5.2          RColorBrewer_1.1-3   
[27] S4Vectors_0.36.2      GlobalOptions_0.1.2  
[29] codetools_0.2-19      shape_1.4.6          
[31] evaluate_0.23         rmarkdown_2.26       
[33] ComplexHeatmap_2.12.1 compiler_4.2.0       
[35] stats4_4.2.0          GetoptLong_1.0.5   
```
