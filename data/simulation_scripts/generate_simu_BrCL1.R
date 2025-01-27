set.seed(0)
today=format(Sys.Date(), format="%y%m%d")

## ----
## Add libs
## ----
source("generic_functions_sim.R")

## ----
## Load data using your own path
## ----
print("-> Loading data...")
T_rna <- readRDS("../references/BrCL1_rna.rds")
T_dnam <- readRDS("../references/BrCL1_dnam.rds")

## ----
## Fix parameters for the simulation
## ----
n_rep = 10
n_samples = 120
varCrit = 10
alpha = c(0.45, 0.1, 0.15, 0.3)
p = .1

## ----
## Generate simulations stored in a folder called 'simulations'
## ----
print(paste0("-> Generate simu..."))
for (i in seq(n_rep)) {
  sim_txt <- ifelse(length(strsplit(as.character(i),"")[[1]])==1,'sim0','sim')
  print(paste0("Simu n°",i))
  data_simu_clean_tot <- generate_simu_tot(alph=alpha,
                                           ref_rna=T_rna,
                                           ref_dnam=T_dnam,
                                           n_samples=n_samples,
                                           varCrit=varCrit)
  Amat <- data_simu_clean_tot$Amat
  Dmat_noise <- add_noise(data_simu_clean_tot,
                         p=p)
  saveRDS(list(D_rna_sim = Dmat_noise$Drna,
               A_ref = Amat),
          file=paste0("../simulations/rna/",today,"_BrCL1_",sim_txt,i,".rds"))
  saveRDS(list(D_dnam_sim = Dmat_noise$Ddnam,
               A_ref = Amat),
          file=paste0("../simulations/dnam/",today,"_BrCL1_",sim_txt,i,".rds"))
}
