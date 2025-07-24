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

prism.states <- function(dat, ref_profiles, nCores = 32) {
  ncores <- nCores - 1
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
  
  ## Deaggregate variable types
  state_labels <- base::colnames(A_state)
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

tpm_norm <- function(dat) {
  TPM = function(counts, lengths) {
    A = intersect(rownames(counts), names(lengths))
    counts = counts[A,]
    lengths = lengths[A]
    rate = counts / lengths
    apply(rate, 2, function(x) 1e6 * x / sum(x))
  }
  human_lengths = readRDS("human_lengths.rds")
  matrix = TPM(counts = as.matrix(dat), lengths = human_lengths)
  rownames(matrix) = toupper(rownames(matrix))
  return(matrix)
}

do_run_sup_deconvolution = function(method, dat, ref_profiles, threads=32) {
  if (method=='InstaPrism') {library(InstaPrism)}
  tictoc::tic(method)
  if (method=="ols") {
    res <- t(granulator::deconvolute(m = tpm_norm(dat), sigMatrix = tpm_norm(ref_profiles), methods = method, use_cores = threads)$
               proportions$
               ols_sig1) #TPM norm
    rownames(res) <- gsub("\\.", "_", rownames(res))
    res <- res[colnames(ref_profiles), colnames(dat)]
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
  }
  else if (method=="fardeep") {
    fardeep = FARDEEP::fardeep(ref_profiles, dat)
    res = t(fardeep$abs.beta)
  }
  else if (method=="fardeepsto") {
    fardeep = FARDEEP::fardeep(ref_profiles, dat)
    res = t(fardeep$relative.beta)
  }
  else if (method=="elasticnet") {
    RESULTS = apply(dat, 2, function(z)
      coef(glmnet::glmnet(x = ref_profiles,
                          y = z,
                          alpha = 0.2,
                          standardize = TRUE,
                          lambda = glmnet::cv.glmnet(as.matrix(ref_profiles), z)$lambda.1se))[1:ncol(ref_profiles) + 1,])
    RESULTS = apply(RESULTS, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    RESULTS = apply(RESULTS, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res = RESULTS
  }
  else if (method=="rlr") { #rlr = robust linear regression
    res <- t(epidish(dat, as.matrix(ref_profiles), method = "RPC")$estF)
  }
  else if (method=="DeconRNASeq") { #NN quadratic programmin
    require(pcaMethods)
    res = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(dat), signatures = as.data.frame(ref_profiles), proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = FALSE, fig = FALSE)$out.all)
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
    colnames(res) = colnames(dat)
  }
  else if (method=="nnls") {
    res <- t(granulator::deconvolute(m = tpm_norm(dat), sigMatrix = tpm_norm(ref_profiles), methods = method, use_cores = threads)$
               proportions$
               nnls_sig1) #TPM norm
    res <- res[gsub("_", ".", colnames(ref_profiles)),]
    rownames(res) <- colnames(ref_profiles)
    res <- apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res <- apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res <- res[, colnames(dat)]
  }
  else if (method=="svr") {
    dat_tpm = tpm_norm(dat)
    res <- t(granulator::deconvolute(m = dat_tpm, sigMatrix = tpm_norm(ref_profiles), methods = 'svr', use_cores = threads)$
               proportions$
               svr_sig1) #TPM norm
    res <- res[gsub("_", ".", colnames(ref_profiles)),]
    rownames(res) <- colnames(ref_profiles)
    res <- res[, colnames(dat)]
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
  }
  else if (method=="CIBERSORTepi") {
    beta.m = dat
    ref.m = as.matrix(ref_profiles)
    res <- t(epidish(beta.m, ref.m, method = "CBS")$estF)
  }
  else if (method == "CIBERSORT") {
    source("scripts/CIBERSORT.R")
    res = CIBERSORT(ref_profiles, dat) 
    res = t(res[,colnames(ref_profiles)])
  }
  else if (method=="WISP") {
    getWeight = function(data, centro, scaling = c("none", "scale", "center")[1], cutoff_gobalFtest = 0.05, Rsquared_cutoff = 0.5, cutoff_ttest_weights = 0.05, sum_LessThanOne = TRUE) {
      g = intersect(rownames(data), rownames(centro))
      data = data[g,]
      centro = as.matrix(centro[g,])
      if (scaling == "scale") {
        data = scale(data, scale = TRUE)
        centro = scale(centro, scale = TRUE)
      }
      if (scaling == "center") {
        data = scale(data, scale = FALSE)
        centro = scale(centro, scale = FALSE)
      }
      nJ = ncol(centro)
      A = centro
      nSubj = ncol(data)
      mixCoef = matrix(0, nSubj, nJ)
      rownames(mixCoef) = colnames(data)
      colnames(mixCoef) = colnames(centro)
      Amat = cbind(rep(-1, nJ), diag(nJ))
      b0vec = c(-1, rep(0, nJ))
      if (sum_LessThanOne == TRUE) {
        meq = 0
      } else {
        meq = 1
      }
      output = data.frame(t(apply(data, 2, function(y) {
        obs = which(!is.na(y))
        Dmat = t(centro[obs,]) %*% centro[obs,]
        diag(Dmat) <- diag(Dmat) + 1e-08
        sc <- norm(Dmat, "2")
        mixCoef = quadprog::solve.QP(Dmat / sc, (t(centro[obs,]) %*%
          y[obs]) / sc, Amat, b0vec, meq = meq)$sol
        B = as.matrix(y)
        coeff = round(mixCoef, 4)
        names(coeff) = paste("weight", colnames(centro), sep = ".")
        vBeta = matrix(coeff)
        dSigmaSq = sum((B - A %*% vBeta)^2) / (nrow(A) - ncol(A))
        dTotalSq = sum((B)^2) / (nrow(A))
        dModelSq = sum((A %*% vBeta)^2) / (ncol(A))
        mVarCovar = try(dSigmaSq * chol2inv(chol(t(A) %*% A)))
        Adjusted.R.squared = round((dTotalSq - dSigmaSq) / dTotalSq, 3)
        Ftest = dModelSq / dSigmaSq
        p.Ftest = stats::pf(q = Ftest, df1 = ncol(A), df2 = nrow(A) - ncol(A), lower.tail = FALSE)
        ng = nrow(centro)
        if (!is.character(mVarCovar)) {
          vStdErr = sqrt(diag(mVarCovar))
          CI.inf = sapply(1:nJ, function(j) { round(coeff[j] - (stats::qt(1 - cutoff_ttest_weights / 2, ng - nJ) * vStdErr[j]), 2) })
          CI.sup = sapply(1:nJ, function(j) { round(coeff[j] + (stats::qt(1 - cutoff_ttest_weights / 2, ng - nJ) * vStdErr[j]), 2) })
          tvalue = sapply(1:nJ, function(j) { coeff[j] / vStdErr[j] })
          p.Ttest = sapply(1:nJ, function(j) { stats::pt(abs(coeff[j] / vStdErr[j]), df = ng - nJ, lower = FALSE) * 2 })
        } else { mVarCovar = NA
          vStdErr = CI.inf = CI.sup = p.Ttest = rep(NA, nJ)
        }
        names(CI.inf) = paste("CI.inf.", colnames(centro), sep = "")
        names(CI.sup) = paste("CI.sup.", colnames(centro), sep = "")
        names(p.Ttest) = paste("Pvalue.", colnames(centro), sep = "")
        coeff.filtered = coeff
        coeff.filtered[p.Ttest > cutoff_ttest_weights] = 0
        names(coeff.filtered) = paste(names(coeff), ".filtered", sep = "")
        dist.Obs.Model = round(sqrt(sum(((centro %*% coeff) - y)^2)), 2)
        c(coeff, dist.Obs.Model = dist.Obs.Model, Ftest.pvalue = p.Ftest, Adjusted.R.squared = Adjusted.R.squared, CI.inf, CI.sup, p.Ttest, coeff.filtered) })))
      output$topWeightedClass = colnames(centro)[apply(output[, 1:nJ], 1, which.max)]
      output$deltaTopWeights = apply(output[, 1:nJ], 1, function(z) abs(diff(sort(z, decreasing = TRUE)[1:2])))
      CI = sapply(1:nJ, function(j) {
        CIinf = sapply(output[, gsub(" ", ".", paste("CI.inf.", colnames(centro)[j], sep = ""))], function(x) max(0, x))
        CIsup = sapply(output[, gsub(" ", ".", paste("CI.sup.", colnames(centro)[j], sep = ""))], function(x) min(1, x))
        paste("[", CIinf, ", ", CIsup, "]", sep = "") })
      colnames(CI) = colnames(centro)
      output$CI = CI
      output$WARNING = as.factor(c("OK", "LIMIT")[(output$Ftest.pvalue > cutoff_gobalFtest | output$Adjusted.R.squared < Rsquared_cutoff) + 1])
      output = output[, -c(grep("CI.inf", colnames(output)), grep("CI.sup", colnames(output)))]
      return(output)
    }
    resW <- getWeight(dat, ref_profiles, scaling = "none")
    res <- t(resW[, grep("filtered", colnames(resW))])
    res <- apply(res, 2, function(x) x / sum(x))
    rownames(res) = colnames(ref_profiles)
  }
  else if (method=="InstaPrism") {
    res <- prism.states(dat, ref_profiles, nCores=threads)
  }
  time_elapsed = tictoc::toc()
  time_elapsed = time_elapsed$toc - time_elapsed$tic
  return(list(res = res, time_elapsed = time_elapsed))
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
  sim_files <- sort(list.files(input_path, pattern = paste0(date, "_", data)))
  # for replicate sim
  sim_file = sim_files[sim]
  sim <- strsplit(strsplit(sim_file, ".rds")[[1]], "_sim")[[1]][[2]]
  data_tot <- readRDS(paste0(input_path, sim_file))
  D <- data_tot[[paste0("D_", omic, "_sim")]]
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
