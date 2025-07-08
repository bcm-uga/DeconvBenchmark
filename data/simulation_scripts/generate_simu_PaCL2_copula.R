set.seed(0)
today=paste0(format(Sys.Date(), format="%y%m%d"),"Copula")

## ----
## Add libs
## ----
source("generic_functions_sim.R")

## ----
## Load data using your own path
## ----
print("-> Loading data...")
T_rna <- readRDS("../references/PaCL2_rna.rds")
T_dnam <- readRDS("../references/PaCL2_dnam.rds")
colnames(T_rna)[2:3] = c("TUM basal","TUM classical")
colnames(T_dnam)[2:3] = c("TUM basal","TUM classical")

empirical_A = readRDS("../invitro/PaMIX_A.rds")
empirical_Drna = readRDS("../invitro/PaMIX_D_rna.rds")
empirical_Ddnam = readRDS("../invitro/PaMIX_D_dnam.rds")
empirical = list(Amat = empirical_A, Drna = empirical_Drna, Ddnam = empirical_Ddnam)

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
                                           ref_dnam=T_dnam,
                                           n_samples=n_samples,
                                           varCrit=varCrit,
                                           dataset_pdac=T)
  Amat <- data_simu_clean_tot$Amat
  rownames(Amat)[2:3] = c("Cancer basal","Cancer classical")
  data_tot<- add_noise_copula(data_simu_clean_tot, empirical, T_rna, T_dnam)
  saveRDS(list(D_rna_sim = data_tot$Drna,
               A_ref = Amat),
          file=paste0("../simulations/rna/",today,"_PaCL2_",sim_txt,i,".rds"))
  saveRDS(list(D_dnam_sim = data_tot$Ddnam,
               A_ref = Amat),
          file=paste0("../simulations/dnam/",today,"_PaCL2_",sim_txt,i,".rds"))
}
