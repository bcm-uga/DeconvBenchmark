set.seed(1)

#####
# Libraries
#####
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
                   "LuCL"=9)

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

do_run_unsup_deconvolution = function(method, dat, Atrue, option=c("Amat","Tmat"), ref_profiles=NULL, dist=F, threads=32) {
  k = nrow(Atrue)
  if (method=="MeDeCom") {library(MeDeCom, quietly=TRUE)}
  if (method=="NMF") {library(NMF)}
  tic(method)
  if (method=="MeDeCom") {
    res <- MeDeCom::runMeDeCom(dat, Ks=k, lambdas=c(0,10^(-4:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=threads, random.seed=1)
    lambda <- res@parameters$lambdas[which.min(res@outputs$`1`$cve)]
    A <- MeDeCom::getProportions(res, K=k, lambda=lambda)
    Tmat <- MeDeCom::getLMCs(res, K=k, lambda=lambda)
    res <- list(A_matrix=A,
                T_matrix=Tmat)
  }
  else if (method=="ICA") {
    res = fastICA::fastICA(X = dat, n.comp = k, maxit = 1000, tol = 1e-09)
    res$names = row.names(dat)
    weighted.list <- deconica::generate_markers(df = res,
                                                n = min(30,nrow(dat)),
                                                return = "gene.ranked")
    ICA_scores_weighted = deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE)
    tmp_dat = t(ICA_scores_weighted)
    colnames(tmp_dat) = colnames(dat)
    res = list(A_matrix=abs(tmp_dat) %*% diag(1/colSums(abs(tmp_dat))),
               T_matrix=res$S)
  }
  else if (method=="NMF") {
    x = dat[rowSums(dat)>0,]
    res <- NMF::nmf(x = x, rank = k, method = "snmf/r", seed = 1)
    A = apply(X = res@fit@H, 2, function(z) {
      z / sum(z) }) #explicit STO constraint
    if (any(rowSums(A) == 0)) {
      A[rowSums(A) == 0, 1] <- 1e-5
    }
    res <- list(A_matrix = A, 
                T_matrix = res@fit@W)
  }
  else if (method=="RefFreeEWAS") {
    res = RefFreeEWAS::RefFreeCellMix(Y = dat, K = k, verbose = F)
    res$Omega = t(apply(res$Omega, 2, function(x)
      ifelse(x < 0, 0, x))) #explicit NN constraints
    res$Omega = apply(res$Omega, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res = list(A_matrix=res$Omega,
               T_matrix=res$Mu)
  }
  else if (method=="EDec") {
    res = t(EDec::run_edec_stage_1(meth_bulk_samples = as.matrix(dat), 
                                   informative_loci = rownames(dat), 
                                   num_cell_types = k)$proportions)
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraints
    res=list(A_matrix=res,
             T_matrix=NULL)
  }
  else if (method=="debCAM") {
    cluster.num = min(5*k, ncol(dat) - 1 )
    if (nrow(dat) < 200) {
      dim.rdc = max(cluster.num, nrow(dat)/10) 
    } else {dim.rdc = 10}
    rCAM <- debCAM::CAM(data = dat,
                K = k,
                cluster.num = cluster.num,
                MG.num.thres = 1,
                lof.thres = 0,
                dim.rdc =  dim.rdc)
    res = list(A_matrix = t(debCAM::Amat(rCAM, k)),
               T_matrix = debCAM::Smat(rCAM, k))
  }
  time_elapsed=toc()
  time_elapsed = time_elapsed$toc - time_elapsed$tic
  #if (option=="Tmat") {
  #  hvg=TOAST::findRefinx(ref_profiles, nmarker = 1e3)
  #  elt1=ref_profiles[hvg,]
  #  elt2=res$T_matrix[hvg,]
  #}
  if (option=="Amat") {
    elt1=t(Atrue)
    elt2=t(res$A_matrix)
  }
  if (dist) {
    row_order <- c(clue::solve_LSAP(dynutils::calculate_distance(t(elt1),t(elt2)), maximum = F))
  }
  else {
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
  }
  res$A_matrix <- res$A_matrix[row_order,]
  rownames(res$A_matrix) <- colnames(elt1)
  return(list(res=res$A_matrix,
              time_elapsed=time_elapsed))
}

SB_deconv_lot_method_sim <- function(lot, block, method, method_class, sim, date, input_path, pred_file, time_file, fs) {
  do_featselec <- ifelse(fs=="none",F,T)
  input_path <- paste0(input_path, block, "/")
  # read files
  T_ref <- as.data.frame(readRDS(paste0(input_path, list.files(input_path, pattern = paste0(date, "_", lot, "_T_", block, "_ref.rds")))))
  sim_files <- sort(list.files(input_path, pattern = paste0(date, "_", lot, "_sim")))
  # for replicate sim
  sim_file = sim_files[sim]
  sim <- strsplit(strsplit(sim_file, ".rds")[[1]], "_sim")[[1]][[2]]
  data <- readRDS(paste0(input_path, sim_file))
  dat <- data[[paste0("D_", block, "_sim")]]
  ref_profiles <- T_ref
  if (!do_featselec) {
    if (block=="met") {
      # restrict size of Dmet to 3e4 features for speed
      hvg_3e4 <- featselec_hvg(dat, n_hvg = 3e4)
      dat <- dat[hvg_3e4,]
      ref_profiles <- ref_profiles[hvg_3e4,]
    }
  }
  if (do_featselec) {
    # proceed to feature selection
    toast_res <- featselec_toast(dat, featselec_K[[lot]])
    dat <- dat[toast_res[[fs]],]
    ref_profiles <- ref_profiles[toast_res[[fs]],]
  }
  # run deconvolution
  Atrue <- data$A_ref
  deconv_res <- do_run_unsup_deconvolution(method, dat, Atrue = Atrue, option = "Amat")
  A_pred <- deconv_res$res
  timing <- deconv_res$time_elapsed
  saveRDS(A_pred, pred_file)
  saveRDS(timing, time_file)
}

#####
# Deconvolution per dataset per method per replicate
#####
args <- commandArgs(trailingOnly = TRUE)
lot = args[1]
block = args[2]
method = args[3]
method_class = args[4]
sim = args[5]
date = args[6]
fs = args[7]
pred_file = args[8]
time_file = args[9]

SB_deconv_lot_method_sim(lot, block, method, method_class,
                         ifelse(sim=="10",10,as.numeric(strsplit(sim,"")[[1]][2])),
                         date, input_path, pred_file, time_file, fs)
