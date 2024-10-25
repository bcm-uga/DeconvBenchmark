This folder contains 
- a folder scripts/ with the R scripts for the deconvolution in each setting (class/omic),
- the Snakefile: this file should be adapted in order to do deconvolution on in vitro/vivo data, or to run the pipeline with new methods,
- the definition file for the apptainer container, along with the directions to create the container from the .def file. The definition file can be modified to test new methods.

In case the user wants to add new datasets, the following parts of the deconvolution scripts need to be modified:
- for the supervised deconvolution with InstaPrism, the function "prism.states" should be modified to define variable types (in our case tumor types) : Sections ##Define variable types and ## Deaggregate variable types,
- the list featselec_K has to include the expected number of cell types k in the new dataset, in the format "New_Data"=k, for the feature selection step,

- in the Snakefile, the parameter DATA_OMIC should be modified as well

In case the user wants to add new methods, he/she should modify:
- the corresponding R script depending if it's a supervised/unsupervised method for RNA/DNAm data,
- the parameter METHOD_OMIC_CLASS in the preamble of the Snakefile,
- the definition file .def of the container if needed.

Tu run the RNA methods that require TPM normalization (OLS, NNLS, SVR), the file with gene lengths "human_lengths.rds" is available upon request from the authors.