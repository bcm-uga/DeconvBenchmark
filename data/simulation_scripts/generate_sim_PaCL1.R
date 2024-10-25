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
T_rna <- readRDS("../references/PaCL1_rna.rds")
T_met <- readRDS("../references/PaCL1_met.rds")

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
                                           ref_met=T_met,
                                           n_samples=n_samples,
                                           varCrit=varCrit,
                                           dataset_pdac=T)
  Amat <- data_simu_clean_tot$Amat
  Dmat_noise <- add_noise(data_simu_clean_tot,
                          p=p)
  saveRDS(list(D_rna_sim = Dmat_noise$Drna,
               A_ref = Amat),
          file=paste0("../simulations/rna/",today,"_PaCL1_",sim_txt,i,".rds"))
  saveRDS(list(D_met_sim = Dmat_noise$Dmet,
               A_ref = Amat),
          file=paste0("../simulations/met/",today,"_PaCL1_",sim_txt,i,".rds"))
}
