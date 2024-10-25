set.seed(0)
today=format(Sys.Date(), format="%y%m%d")

## ----
## Add libs
## ----
lot="dPANCREAS"
source(paste0("../../src/0simu_dirichlet_functions_",lot,".R"))

## ----
## load data
## ----
print("-> Loading data...")
T_rna <- readRDS("../../../datashare/DECONbench/PANCREAS/results/T_rna_raw.rds")
T_met <- readRDS("../../../datashare/DECONbench/PANCREAS/results/T_met_raw.rds")
saveRDS(T_rna, file=paste0("simulations/rna/",today,"_",lot,"_T_rna_ref.rds"))
saveRDS(T_met, file=paste0("simulations/met/",today,"_",lot,"_T_met_ref.rds"))

## ----
## generate simulated set
## ----
n_samples=120
params = data.frame(varCrit=10,sd_rna=1,sd_met=3)
print(paste0("-> Generate training set..."))
n_rep=10
for (i in seq(n_rep)) {
  sim_txt <- ifelse(length(strsplit(as.character(i),"")[[1]])==1,'sim0','sim')
  if (!file.exists(paste0("simulations/met/",today,"_",lot,"_",sim_txt,i,".rds")) &
      !file.exists(paste0("simulations/rna/",today,"_",lot,"_",sim_txt,i,".rds")))
  print(paste0("Simu n°",i))
  data_simu_clean0 <- generate_simu_set(T_rna,T_met, n_samples=n_samples, varCrit=params$varCrit)
  A <- data_simu_clean0[[3]]
  data_simu_clean <- data_simu_clean0[1:2]
  data_simu <- generate_simu_noise(data_simu_clean$D_rna,data_simu_clean$D_met,
                                   p=0.1,sd_rna=params$sd_rna,sd_met=params$sd_met) # data raw
  saveRDS(list(D_rna_sim=data_simu[[1]],
               A_ref=A), file=paste0("simulations/rna/",today,"_",lot,"_",sim_txt,i,".rds"))
  saveRDS(list(D_met_sim=data_simu[[2]],
               A_ref=A), file=paste0("simulations/met/",today,"_",lot,"_",sim_txt,i,".rds"))
}
