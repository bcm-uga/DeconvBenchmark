# A robust workflow to benchmark deconvolution of multi-omic data

This repository contains all the files needed to perform deconvolution on example datasets, along with the ranking and the figures of the paper.
Our pipeline can be easily extended to include and evaluate novel methods, as well as other datasets.

## How to make in silico data: folder [data](data/)

Shortly, we made in silico data using reference profiles of pure cell types from different tissues convoluted with proportions generated based on a Dirichlet distribution. The scripts are in [simulation_scripts](data/simulation_scripts/), and the reference profiles can be downloaded from Zenodo DOIDOIDOI.
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

This is also where we stored the *in vitro* and *in vivo* datasets, with the syntax data/invitro/DATA_A.rds for the reference matrix and data/invitro/DATA_D_OMIC.rds for the methylation/expression matrix/. Please refer to the table [here](data/README.md) for the instructions on where to download *in vitro* and *in vivo* data.

## How to run the deconvolution methods: folder [deconvolution](deconvolution/)

This pipeline uses an Apptainer container and Snakemake. Instructions on how to use these tools can be found [here](deconvolution/README.md) and [here](deconvolution/container.md). The Apptainer container can be downloaded from Zenodo DOIDOIDOI.

Briefly, this folder contains the scripts to perform the deconvolution pipeline. There is one script per setting (class of the method/omic type). Those scripts can be modified to include new methods and/or new datasets (cf [README](deconvolution/README.md)).

Snakemake will run all methods for all omics. The [Snakefile](deconvolution/Snakefile) is self-explanatory and can be modified to include new methods/datasets. In general, you can refer to the [README](deconvolution/README.md) to know how to test new methods/datasets.

Results of the deconvolution, *i.e.* estimation of the proportion matrix along with elapsed time will be stored in deconvolution/results/prediction/OMIC/CLASS/ for the proportion matrix and deconvolution/results/timing/OMIC/CLASS/ for the time elapsed with the syntax 240101_DATA1_Apred_FS_METHOD_sim01.rds / 240101_DATA1_timing_FS_METHOD_sim01.rds (FS being the feature selection strategy and METHOD the deconvolution algorithm).

## How to do the ranking and reproduce the figures of the paper: folder [ranking_figures](ranking_figures/)

