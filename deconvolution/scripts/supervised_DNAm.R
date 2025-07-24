set.seed(1)

#####
# Libraries
#####
library(EpiDISH)
library(tictoc)

#####
# Set parameters, change with your own path
#####
input_path <- "../data/simulations/"
# true number of cell types in each dataset for the TOAST feature selection step
featselec_K = list("BrCL1"=4,
                   "BrCL2"=6,
                   "PaCL1"=5,
                   "PaCL2"=9,
                   "BlCL"=6,
                   "LuCL"=9,
                   "PaPB"=7)

#####
# Functions
#####
featselec_hvg <- function(dat, n_hvg) {
  hvg <- TOAST::findRefinx(dat, nmarker = n_hvg)
  return(hvg)
}

featselec_toast <- function(dat, k) {
  toast <- TOAST::csDeconv(dat, k, TotalIter = 10, FUN = function (dat.arg, k.arg) {
    res <- fastICA::fastICA(X = dat.arg, n.comp = k.arg, maxit = 1000, tol = 1e-09)
    res$names <- row.names(dat.arg)
    weighted.list <- deconica::generate_markers(df = res, n = 30, return = "gene.ranked")
    res$A_weighted <- t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
    colnames(res$A_weighted) <- colnames(dat.arg)
    A <- abs(res$A_weighted) %*% diag(1 / colSums(abs(res$A_weighted)))
    return(t(A))
  })$updatedInx
  hvf <- TOAST::findRefinx(dat, nmarker = length(toast))
  return(list(toast = toast, hvf = hvf))
}

prism.states <- function(dat, ref_profiles, nCores = threads) {
  ncores <- nCores - 1
  dat <- prism.met.to.count(dat) # transform DNAm data to pseudo counts
  ref_profiles <- prism.met.to.count(ref_profiles) # transform DNAm data to pseudo counts
  # define types and states
  state_labels <- colnames(ref_profiles)
  type_labels <- state_labels
  
  ## Define variable types (tumor types in our case)
  ## BrCL1: types and states are equals, tumoral type label is "tumor"
  ## PaCL1: 2 tumoral states, 1 tumoral type
  type_labels[grepl("TUM_", type_labels)] <- "tumor"
  ## PaCL2 and PaPB: 2 tumoral states, 1 tumoral type
  type_labels[grepl("Cancer", type_labels)] <- "tumor"
  ## LuCL: 1 tumoral state and type
  type_labels[grepl("A549", type_labels)] <- "tumor"
  ## BrCL2: 3 tumoral states, 1 tumoral type
  type_labels[grepl("BT474", type_labels)] <- "tumor"
  type_labels[grepl("MCF7", type_labels)] <- "tumor"
  type_labels[grepl("T47D", type_labels)] <- "tumor"
  
  
  prism <- BayesPrism::new.prism(
      reference = base::t(ref_profiles),
      mixture = base::t(dat),
      input.type = "GEP",
      cell.state.labels = state_labels,
      cell.type.labels = type_labels,
      key = "tumor") # create prism obj
  
  res <- InstaPrism(prismObj = prism, input_type = "prism",
                    n.core = ncores) # run deconv
  A_state <- t(res@Post.ini.cs@theta) # get state props before update
  A_type <- t(res@Post.ini.ct@theta) # get type props after update
  # get results
  ## we need to "deaggregate" tumoral type to tumoral states
  ### first get tumoral states proportions within the tumoral type
  state_labels <- base::colnames(A_state)
  
  ## Deaggregate variable types
  tumoral_states_mask <- state_labels == "tumor" |
    grepl("TUM_", state_labels) |
    grepl("Cancer", state_labels) |
    grepl("A549", state_labels) |
    grepl("BT474", state_labels) |
    grepl("MCF7", state_labels) |
    grepl("T47D", state_labels)
  
  tumoral_states_labels <- state_labels[tumoral_states_mask] # use of labels rather than mask to make sure order in the dimensions name is not important
  not_tumoral_states_labels <- state_labels[!tumoral_states_mask]
  tumoral_states_prop <- A_state[, tumoral_states_labels]
  if (!is.array(tumoral_states_prop)) {
    tumoral_states_prop <- matrix(tumoral_states_prop, ncol = 1,
                                  dimnames = list(base::rownames(A_state), tumoral_states_labels))
  }
  tumoral_states_prop <- tumoral_states_prop / apply(tumoral_states_prop, 1, sum)
  #### then dispatch tumoral type final proportion to states accordingly
  type_labels <- base::colnames(A_type)
  tumoral_type_mask <- type_labels == "tumor"
  tumoral_type_labels <- type_labels[tumoral_type_mask]
  not_tumoral_type_labels <- type_labels[!tumoral_type_mask]
  A_matrix <- array(dim = dim(A_state), dimnames = dimnames(A_state))
  A_matrix[, tumoral_states_labels] <- A_type[, tumoral_type_labels] * tumoral_states_prop
  ### other states remain unchanged
  A_matrix[, not_tumoral_states_labels] <- A_type[, not_tumoral_type_labels]
  return(t(A_matrix[, colnames(ref_profiles)]))
}

prism.met.to.count <- function(dat, factor = 1000) {
  return(round(dat * factor))
}

do_run_sup_deconvolution <- function(method, dat, ref_profiles, threads=32) {
  if (method == 'InstaPrism') {library(InstaPrism)}
  tictoc::tic(method)
  if (method=="rlr") {
    beta.m = dat
    ref.m = as.matrix(ref_profiles)
    rownames(ref.m) = rownames(beta.m)
    res <- t(epidish(beta.m, ref.m, method = "RPC")$estF)
  }
  else if (method=="nnls") {
    res <- t(granulator::deconvolute(m = dat, sigMatrix = as.matrix(ref_profiles), methods = method, use_cores = threads)$
               proportions$
               nnls_sig1)
    res <- res[gsub("_", ".", colnames(ref_profiles)),]
    rownames(res) <- colnames(ref_profiles)
    res <- apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res <- apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res <- res[, colnames(dat)]
  }
  else if (method=="CIBERSORTepidish") {
    beta.m = dat
    ref.m = as.matrix(ref_profiles)
    res <- t(epidish(beta.m, ref.m, method = "CBS")$estF)
  }
  else if (method == "CIBERSORT") {
    source("scripts/CIBERSORT.R")
    res = CIBERSORT(ref_profiles, dat) 
    res = t(res[,colnames(ref_profiles)])
  }
  else if (method=="epidishCP") {
    res <- t(epidish(beta.m = dat, ref.m = as.matrix(ref_profiles), method = "CP")$estF)
    res[res<0]<-0
  }
  else if (method=="InstaPrism") {
    res <- prism.states(dat, ref_profiles, nCores = threads)
  }
  time_elapsed = tictoc::toc()
  time_elapsed = time_elapsed$toc - time_elapsed$tic
  return(list(res=res,time_elapsed=time_elapsed))
}

SB_deconv_data_method_sim <- function(data, omic, method, method_class, sim, date, input_path, pred_file, time_file, fs) {
  do_featselec <- ifelse(fs=="none",F,T)
  input_path_ref <- paste0(input_path,"../references/")
  input_path <- paste0(input_path, omic, "/")
  # read files
  list_files = list.files(input_path_ref, pattern = paste0(data), full.names = T)
  if (length(grep("_sc",list_files))>0) {list_files = list_files[-grep("_sc",list_files)]}
  if (length(list_files)>1) {list_files = grep(omic,list_files,value = T)}
  ref_profiles <- as.data.frame(readRDS(list_files))
  sim_files <- sort(list.files(input_path, pattern = paste0(date, "_", data, "_sim")))
  # for replicate sim
  sim_file = sim_files[sim]
  sim <- strsplit(strsplit(sim_file, ".rds")[[1]], "_sim")[[1]][[2]]
  data_tot <- readRDS(paste0(input_path, sim_file))
  D <- data_tot[[paste0("D_", omic, "_sim")]]
  if (!do_featselec) {
    if (omic=="dnam") {
      # restrict size of Ddnam and ref matrix to 3e4 features for speed
      hvg_3e4 <- featselec_hvg(D, n_hvg = 3e4)
      D <- D[hvg_3e4,]
      ref_profiles <- ref_profiles[hvg_3e4,]
    }
  }
  if (do_featselec) {
    # proceed to feature selection
    toast_res <- featselec_toast(D, featselec_K[[data]])
    D <- D[toast_res[[fs]],]
    ref_profiles <- ref_profiles[toast_res[[fs]],]
  }
  # run deconvolution
  ref_profiles <- ref_profiles[, sort(colnames(ref_profiles))]
  deconv_res <- do_run_sup_deconvolution(method, D, ref_profiles)
  A_pred <- deconv_res$res
  timing <- deconv_res$time_elapsed
  saveRDS(A_pred, pred_file)
  saveRDS(timing, time_file)
}
                   
#####
# Deconvolution per dataset per method per replicate
#####
args <- commandArgs(trailingOnly = TRUE)
data = args[1]
omic = args[2]
method = args[3]
method_class = args[4]
sim = args[5]
date = args[6]
fs = args[7]
pred_file = args[8]
time_file = args[9]

SB_deconv_data_method_sim(data, omic, method, method_class,
                         ifelse(sim=="10",10,as.numeric(strsplit(sim,"")[[1]][2])),
                         date, input_path, pred_file, time_file, fs)
