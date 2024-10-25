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
T_rna <- readRDS("../references/PaCL2_rna.rds")
T_met <- readRDS("../references/PaCL2_met.rds")
colnames(T_rna)[2:3] = c("TUM basal","TUM classical")
colnames(T_met)[2:3] = c("TUM basal","TUM classical")

## ----
## Fix parameters for the simulation
## ----
n_rep = 10
n_samples = 120
varCrit = 10
alpha = c(0.01, 0.04, 0.03, 0.15, 0.46, 0.01, 0.01, 0.29)
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
  rownames(Amat)[2:3] = c("Cancer basal","Cancer classical")
  Dmat_noise <- add_noise(data_simu_clean_tot,
                          p=p)
  saveRDS(list(D_rna_sim = Dmat_noise$Drna,
               A_ref = Amat),
          file=paste0("../simulations/rna/",today,"_PaCL2_",sim_txt,i,".rds"))
  saveRDS(list(D_met_sim = Dmat_noise$Dmet,
               A_ref = Amat),
          file=paste0("../simulations/met/",today,"_PaCL2_",sim_txt,i,".rds"))
}
