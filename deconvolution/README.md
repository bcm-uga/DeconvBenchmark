This folder contains 
- a folder scripts/ with the R scripts for the deconvolution in each setting (class/omic),
- the Snakefile: this file should be adapted in order to do deconvolution on in vitro/vivo data, or to run the pipeline with new methods,
- the definition file for the apptainer container, along with the directions to create the container from the .def file (container.md). The definition file can be updated to test new methods by adding necessary packages.

## Quick guide for snakemake install with mamba
```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n YOUR_ENV 'snakemake==7.22.0'
conda activate YOUR_ENV
```

## Run R scripts from a Snakefile file
- The parameter 'input_path' in the R scripts needs to be changed to the path where you stored the data to be deconvoluted.
- For the snakefile, the path YOUR/PROJECT/ROOT needs to be changed to the root path of your project in the shell command of the deconvoution rules.
- Run the Snakefile (here with the OAR job manager)
```
snakemake --cluster "oarsub --project YOUR_PROJECT -l /nodes=1/core=1,walltime=24:00:00" --latency-wait 60 --cores 1 --jobs 50
```

## Adding new datasets
In case the user wants to add new datasets, the following parts of the deconvolution scripts need to be modified:
- for the supervised deconvolution with InstaPrism, the function "prism.states" should be modified to include the variable cell types (in our case tumor types) of the new dataset: sections ##Define variable types and ## Deaggregate variable types,
- the list featselec_K has to include the expected number of cell types k in the new dataset, in the format "New_Data"=k, for the feature selection step,
- in the Snakefile, the parameter DATA_OMIC should be modified as well

## Adding new methods
In case the user wants to add new methods, he/she should modify:
- the corresponding R script depending on if it's a supervised/unsupervised method for RNA/DNAm data,
- the parameter METHOD_OMIC_CLASS in the preamble of the Snakefile,
- the definition file .def of the container if needed to add required packages to the container.

Tu run the RNA methods that require TPM normalization (OLS, NNLS, SVR), the file with gene lengths "human_lengths.rds" is available upon request from the authors.