## Contents
This folder contains:
- a folder [scripts](scripts) with the R scripts for the deconvolution in each setting (class/omic),
- the [Snakefile](Snakefile): this file should be adapted in order to do deconvolution on in vitro/vivo data or new datasets, or to run the pipeline with new methods (see below),
- the definition file for the apptainer container [container2.def](container2.def), along with the [directions](container.md) to create the container file .sif from the definition file. The definition file can be updated to test new methods by adding necessary packages.

## Quick guide to install snakemake with mamba
```shell
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n YOUR_ENV 'snakemake==7.22.0'
conda activate YOUR_ENV
```

## Launch the Snakemake workflow
- The parameter 'input_path' in the R scripts needs to be changed to the path where you stored the dataset to be deconvoluted.
- In [Snakefile](Snakefile), the path YOUR/PROJECT/ROOT needs to be changed to the root path of your project in the shell command of the deconvoution rules, as well as the parameter DATE to match the date of the day you did your simulations.
- Run the Snakefile
```shell
snakemake --latency-wait 60 --cores 1 --jobs 50
```

## Adding new datasets
In case the user wants to add new datasets, the following parts of the deconvolution scripts need to be modified:
- for the supervised deconvolution with InstaPrism, the function "prism.states" in [scripts](scripts/) should be modified to include the variable cell types (in our case tumor types) of the new dataset: sections ##Define variable types and ##Deaggregate variable types,
- the list featselec_K in the [scripts](scripts/) has to include the expected number of cell types k in the new dataset, in the format "New_Data"=k, for the feature selection step,
- in the [Snakefile](Snakefile), the parameter DATA_OMIC should be modified as well to include the new dataset.

## Adding new methods
In case the user wants to add new methods, he/she should modify:
- the corresponding R script depending on if it's a supervised/unsupervised method for RNA/DNAm data,
- the parameter METHOD_OMIC_CLASS in the preamble of the [Snakefile](Snakefile),
- the definition file [container2.def](container2.def) of the container if needed to add required packages to the container.

To run the RNA methods that require TPM normalization (OLS, NNLS, SVR), the file with gene lengths "human_lengths.rds" is available upon request from the authors.
