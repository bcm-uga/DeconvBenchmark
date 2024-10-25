set.seed(0)
today=format(Sys.Date(), format="%y%m%d")

## ----
## Add libs
## ----
source("generic_functions.R")

## ----
## Load data using your own path
## ----
print("-> Loading data...")
T_rna <- readRDS("../references/BlCL.rds")

## ----
## Fix parameters for the simulation
## ----
n_rep=10
n_samples=120
params = list(varCrit=10,sd_rna=1)
alpha

## ----
## Generate simulations stored in a folder called 'simulations'
## ----
print(paste0("-> Generate simu..."))
for (i in seq(n_rep)) {
  sim_txt <- ifelse(length(strsplit(as.character(i),"")[[1]])==1,'sim0','sim')
  print(paste0("Simu n°",i))
  data_simu_clean_tot <- generate_simu_set(T_rna,
                                           n_samples=n_samples,
                                           varCrit=params$varCrit)
  Amat <- data_simu_clean_tot[[2]]
  Dmat <- data_simu_clean_tot[[1]]
  Dmat_noise <- generate_simu_noise(Dmat,
                                    p=0.1,
                                    sd_rna=params$sd_rna)
  saveRDS(list(D_rna_sim = Dmat_noise,
               A_ref = Amat),
          file=paste0("../simulations/rna/",today,"_BlCL_",sim_txt,i,".rds"))
}

