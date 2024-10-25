set.seed(0)
today <- format(Sys.Date(), format = "%y%m%d")

## ----
## Add libs
## ----
source(paste0("../../src/0simu_load_data_lot.R"))
source(paste0("../../src/0simu_dirichlet_functions_lot.R"))

lot <- "lot1"

## ----
## load data
## ----
print("-> Loading data...")
data_gt <- load_data_lot(lot)
data <- data_gt$data
gt <- data_gt$gt
n_cell_types <- nrow(gt)
ref_row_order <- unname(apply(gt[, (ncol(gt) - n_cell_types + 1):ncol(gt)],
                              2,
                              function(x) rownames(gt)[x == 1]))
T_rna <- data$rna[, (ncol(data$rna) - n_cell_types + 1):ncol(data$rna)]
colnames(T_rna) <- ref_row_order
T_met <- data$met[, (ncol(data$met) - n_cell_types + 1):ncol(data$met)]
colnames(T_met) <- ref_row_order
celltypes <- rownames(gt)
rm(gt)
saveRDS(T_rna[, sort(colnames(T_rna))], file = paste0("simulations/rna/", today, "_", lot, "_T_rna_ref.rds"))
saveRDS(T_met[, sort(colnames(T_met))], file = paste0("simulations/met/", today, "_", lot, "_T_met_ref.rds"))

## ----
## generate simulated set
## ----
n_samples <- 120
print(paste0("-> Generate training set..."))
n_rep <- 10
for (i in seq(n_rep)) {
  sim_txt <- ifelse(length(strsplit(as.character(i), "")[[1]]) == 1, 'sim0', 'sim')
  if (!file.exists(paste0("simulations/met/", today, "_", lot, "_", sim_txt, i, ".rds")) &
    !file.exists(paste0("simulations/rna/", today, "_", lot, "_", sim_txt, i, ".rds")))
    print(paste0("Simu n°", i))
  data_simu_clean0 <- generate_simu_set(data, celltypes, n_samples = n_samples, varCrit = 10)
  prop_simu <- data_simu_clean0[[3]]
  data_simu_clean <- data_simu_clean0[1:2]
  data_simu <- generate_simu_noise(data_simu_clean[[1]], data_simu_clean[[2]],
                                   p = 0.1) # data raw
  saveRDS(list(D_rna_sim = data_simu[[1]],
               A_ref = prop_simu), file = paste0("simulations/rna/", today, "_", lot, "_", sim_txt, i, ".rds"))
  saveRDS(list(D_met_sim = data_simu[[2]],
               A_ref = prop_simu), file = paste0("simulations/met/", today, "_", lot, "_", sim_txt, i, ".rds"))
}
