This folder contains 
- a folder scripts/ with the R scripts for each setting (class/omic) for the deconvolution pipeline,
- the Snakefile: this file should be adapted in order to do deconvolution on in vitro/vivo data,
- the definition file for the apptainer container, along with the directions to create the container from the .def file.

In case the user wants to add new datasets, the following parts of the deconvolution scripts need to be modified:
- for the supervised deconvolution with InstaPrism, the function "prism.states" should be modified to define variable types (in our case tumor types) : Sections ##Define variable types and ## Deaggregate variable types,
- the list featselec_K has to include the expected number of cell types k in the new dataset, in the format "New_Data"=k, for the feature selection step

For the RNA methods that require TPM normalization (OLS, NNLS, SVR), the file "human_lengths.rds" is available upon request from the authors.