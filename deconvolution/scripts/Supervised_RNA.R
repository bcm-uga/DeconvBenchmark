set.seed(1)

# Functions
library(EpiDISH)
library(tictoc)

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

prism.states <- function(dat, ref_profiles, methy = FALSE, nCores = 32) {
  ncores <- nCores - 1
  if (methy) {
    # transform met data to "pseudo"-count data
    dat <- prism.met.to.count(dat)
  }
  # define types and states
  state_labels <- colnames(ref_profiles)
  type_labels <- state_labels
  ## dBREAST: types and states are equals, tumoral type label is "tumor"
  ## dPANCREAS: 2 tumoral states, 1 tumoral type
  type_labels[grepl("TUM_", type_labels)] <- "tumor"
  ## lot1: 2 tumoral states, 1 tumoral type
  type_labels[grepl("Cancer ", type_labels)] <- "tumor"
  ## He: 1 tumoral state and type
  type_labels[grepl("A549", type_labels)] <- "tumor"
  ## Cobos: 3 tumoral states, 1 tumoral type
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
  tumoral_states_mask <- state_labels == "tumor" |
    grepl("TUM_", state_labels) |
    grepl("Cancer ", state_labels) |
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

tpm_norm <- function(dat) {
  TPM = function(counts, lengths) {
    A = intersect(rownames(counts), names(lengths))
    counts = counts[A,]
    lengths = lengths[A]
    rate = counts / lengths
    apply(rate, 2, function(x) 1e6 * x / sum(x))
  }
  human_lengths = readRDS("/bettik/PROJECTS/pr-epimed/amblaeli/projects/acacia/src/TPM/human_lengths.rds")
  matrix = TPM(counts = as.matrix(dat), lengths = human_lengths)
  rownames(matrix) = toupper(rownames(matrix))
  return(matrix)
}

do_run_sup_deconvolution = function(method, dat, ref_profiles, threads=32) {
  if (method == 'InstaPrism') {library(InstaPrism)}
  tic(method)
  if (method == "ols") {
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
  else if (method == "fardeep") {
    fardeep = FARDEEP::fardeep(ref_profiles, dat)
    res = t(fardeep$abs.beta)
  }
  else if (method == "fardeepsto") {
    fardeep = FARDEEP::fardeep(ref_profiles, dat)
    res = t(fardeep$relative.beta)
  }
  else if (method == "DeconFeature2") {
    DeconFeature_fun <- function(x, y, tol = 1e-08, max.iter = 1000, learning_rate = 0.6) {
      m <- nrow(x)
      p <- ncol(x)
      mx <- colMeans(x)
      cx <- as.matrix(x-tcrossprod(rep(1,m),mx))
      cy <- y-mean(y)
      beta <- rep(1/p, p)  # uniform
      fit <- cx%*%beta
      e <- (cy-fit)[,1]
      objective <- sum(e^2)
      iter <- 0
      delta <- 1
      while((delta>tol) & (iter<max.iter)) {
        iter <- iter+1
        current_b <- beta
        m_b <- matrix(0,nrow=p,ncol=p)
        v_obj <- rep(0,p)
        for (j in 1:p) {
          b <- current_b
          fit <- cx%*%b
          e <- (cy-fit)[,1]
          grad_j <- -2*crossprod(cx[,j],e) * b[j]*(1-b[j])
          b[j] = b[j] - learning_rate*grad_j
          b[-j] <- b[-j]*(1-b[j])/sum(b[-j]) # sum-to-one constraint
          b3 <- b
          fit <- cx%*%b
          e <- (cy-fit)[,1]
          obj3 <- sum(e^2)
          b <- current_b
          b4 <- b
          fit <- cx%*%b
          e <- (cy-fit)[,1]
          obj4 <- sum(e^2)
          obj <- c(obj3,obj4)
          b <- cbind(b3,b4)
          which_min <- which.min(obj)
          v_obj[j] <- obj[which_min]
          m_b[,j] <- b[,which_min]
        }
        j_opt <- which.min(v_obj)
        beta <- m_b[,j_opt]
        objective <- min(v_obj)
        new_b <- beta
        delta <- mean(abs(current_b-new_b))
        current_b <- new_b
      }
      fit <- cx%*%beta
      e <- (cy-fit)[,1]
      obj <- sum(e^2)
      trace_beta <- matrix(beta,ncol=1)
      iter <- 0
      delta <- 1
      nonzero <- which(abs(beta)>1e-08)
      l <- log(beta[nonzero]/(1-beta[nonzero]))
      Sx <- crossprod(cx[,nonzero])/(m-1)
      iSx <- solve(Sx)
      while((delta>tol)&(iter<max.iter)&(length(nonzero)>1)) {
        iter <- iter+1
        fit <- cx%*%beta
        e <- (cy-fit)[,1]
        objective <- sum(e^2)
        obj <- c(obj,objective)
        trace_beta <- cbind(trace_beta,beta)
        score <- -2*crossprod(cx[,nonzero],e)
        new_l <- l-(diag(1/(beta[nonzero]*(1-beta[nonzero])),
                        nrow=length(nonzero),
                        ncol=length(nonzero))%*%iSx%*%score/(2*(m-1)))[,1]
        new_beta <- rep(0,times=p)
        new_beta[nonzero] <- 1/(1+exp(-new_l))
        new_beta <- new_beta/sum(new_beta)
        nonzero <- which(abs(new_beta)>1e-08)
        new_l <- log(new_beta[nonzero]/(1-new_beta[nonzero]))
        Sx <- crossprod(cx[,nonzero])/(m-1)
        iSx <- solve(Sx)
        delta <- mean(abs(beta-new_beta), na.rm=T)
        beta <- new_beta
        l <- new_l
      }
      beta <- trace_beta[,which.min(obj),drop=TRUE]
      converged <- TRUE
      if (iter>=max.iter) {
        converged <- FALSE
        warning("Maximum number of iterations reached")
      }
      l <- log(beta/(1-beta))
      fit <- cx%*%beta
      epsilon <- (cy-fit)[,1]
      rmse <- sqrt(mean(epsilon^2))
      return(list(
        l = l,
        beta = beta,
        converged = converged,
        fit = fit,
        rmse = rmse
      ))
    }
    DeconF_test = lapply(seq(ncol(dat)), function(sample) {
      DeconFeature_fun(ref_profiles, dat[,sample])})
    res = matrix(sapply(seq(ncol(dat)), function(sample) {DeconF_test[[sample]]$beta}), 
                 nrow = ncol(ref_profiles), 
                 dimnames = list(colnames(ref_profiles), paste0("train",seq(ncol(dat)))))
  }
  else if (method == "elasticnet") {
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
  else if (method == "rlr") { #rlr = robust linear regression
    res <- t(epidish(dat, as.matrix(ref_profiles), method = "RPC")$estF)
  }
  else if (method == "DeconRNASeq") { #NN quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
    require(pcaMethods)
    res = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(dat), signatures = as.data.frame(ref_profiles), proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = FALSE, fig = FALSE)$out.all)
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
    colnames(res) = colnames(dat)
  }
  else if (method == "nnls") {
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
  else if (method == "svr3") {
    dat_tpm = tpm_norm(dat)
    hvg3_tpm = featselec_hvg(dat_tpm, n_hvg =  min(1e3,nrow(dat_tpm)))
    res <- t(granulator::deconvolute(m = dat_tpm[hvg3_tpm,], sigMatrix = tpm_norm(ref_profiles)[hvg3_tpm,], methods = 'svr', use_cores = threads)$
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
  else if (method == "svr") {
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
  else if (method == "CIBERSORT3") {
    hvg3 = featselec_hvg(dat, n_hvg = min(1e3,nrow(dat)))
    beta.m = dat[hvg3,]
    ref.m = as.matrix(ref_profiles[hvg3,])
    res <- t(epidish(beta.m, ref.m, method = "CBS")$estF)
  }
  else if (method == "CIBERSORT") {
    beta.m = dat
    ref.m = as.matrix(ref_profiles)
    res <- t(epidish(beta.m, ref.m, method = "CBS")$estF)
  }
  else if (method == "WISP") {
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
  else if (method == "InstaPrism") {
    res <- prism.states(dat, ref_profiles, nCores=threads)
  }
  time_elapsed = toc()
  time_elapsed = time_elapsed$toc - time_elapsed$tic
  return(list(res = res, time_elapsed = time_elapsed))
}

SB_deconv_lot_method_sim <- function(lot, block, method, method_class, sim, date, input_path, pred_file, time_file, fs) {
  do_featselec <- ifelse(fs=="none",F,T)
  input_path <- paste0(input_path, block, "/")
  # read files
  T_ref <- as.data.frame(readRDS(paste0(input_path, list.files(input_path, pattern = paste0(date, "_", lot, "_T_", block, "_ref.rds")))))
  sim_files <- sort(list.files(input_path, pattern = paste0(date, "_", lot, "_sim")))
  # for simulation sim
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
  ref_profiles <- ref_profiles[, sort(colnames(ref_profiles))]
  deconv_res <- do_run_sup_deconvolution(method, dat, ref_profiles)
  A_pred <- deconv_res$res
  timing <- deconv_res$time_elapsed
  saveRDS(A_pred, pred_file)
  saveRDS(timing, time_file)
}

## ----
## Set parameters
## ----
input_path <- "/bettik/PROJECTS/pr-epimed/amblaeli/projects/acacia_2final/results/0simu/simulations/"
featselec_K = list("dBREAST"=4,
                   "dPANCREAS"=5,
                   "lot1"=9,
                   "Cobos"=6,
                   "Hoek"=6,
                   "He"=9)
                   
## ----
## Deconvolution per lot per method per sim
## ----
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
