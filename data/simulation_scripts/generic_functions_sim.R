generate_proportions <- function(n_samples, celltypes, alph, varCrit, dataset_pdac) {
  # Compute proportions of various compartments
  if (! dataset_pdac) {
    comp_mix <-t(gtools::rdirichlet(n = n_samples, alpha=alph*varCrit))
    rownames(comp_mix) = celltypes
  }
  else {
    num.cell.types = length(celltypes)
    # Compute proportions of all compartments
    global_mix <-t(gtools::rdirichlet(n = n_samples, alpha=alph*varCrit))
    # Compute relative proportions of tumor types with varying proportions
    #tumor_mix <- t(gtools::rdirichlet(n = n_samples, alpha = c(0.5, 0.5)))
    tumor_mix <- cbind(t(gtools::rdirichlet(n = 1*n_samples / 6, alpha = c(0.5, 0.5))),
                       t(gtools::rdirichlet(n = 2*n_samples / 6, alpha = c(0.25, 0.75))),
                       t(gtools::rdirichlet(n = 1*n_samples / 6, alpha = c(0.75, 0.25))),
                       t(gtools::rdirichlet(n = 2*n_samples / 6, alpha = c(0.1, 0.9))))
    # Merge total matrix proportion
    comp_mix <- matrix(nrow = num.cell.types, ncol = n_samples)
    rownames(comp_mix) = c(celltypes[-grep("TUM",celltypes)],
                          celltypes[grep("TUM",celltypes)])
    comp_mix[1:(num.cell.types - 2), ] <- global_mix[1:(num.cell.types - 2),]
    comp_mix[num.cell.types - 1, ] <- tumor_mix[1, ] * global_mix[num.cell.types - 1,]
    comp_mix[num.cell.types,] <- tumor_mix[2, ] * global_mix[num.cell.types - 1,]
    comp_mix <- comp_mix[celltypes,]
  }
  colnames(comp_mix) = paste("train", 1:n_samples, sep = "")
  return(comp_mix)
}

generate_simu_tot <- function(alph, ref_rna=NULL, ref_dnam=NULL, n_samples=120, varCrit=10, dataset_pdac=F) {
  if (!is.null(ref_rna)) {
    celltypes <- colnames(ref_rna)
  } else {celltypes <- colnames(ref_dnam)}
  Amat <- generate_proportions(n_samples, celltypes, alph, varCrit, dataset_pdac) #prop x sample
  # convolution
  if (!is.null(ref_rna)) {
    Drna <- as.matrix(ref_rna[,celltypes]) %*% Amat[celltypes,]
  }
  if (!is.null(ref_dnam)) {
    Ddnam <- as.matrix(ref_dnam[,celltypes]) %*% Amat[celltypes,]
  }
  result = list()
  if (!is.null(ref_rna)) {
    result$Drna <- Drna
  }
  if (!is.null(ref_dnam)) {
    result$Ddnam <- Ddnam
  }
  result$Amat <- Amat[celltypes,]
  return(result)
}

add_noise <- function(result, p, sd_rna=1, sd_dnam=3) {
  omic = grep("D",names(result),value=T)
  result_noise = list()
  if ("Drna" %in% omic) {
    Drna_noise <- add_noise_nb(result$Drna, p, sd_rna)
    rownames(Drna_noise)=rownames(result$Drna)
    colnames(Drna_noise)=colnames(result$Drna)
    result_noise$Drna <- Drna_noise
  }
  if ("Ddnam" %in% omic) {
    beta_val <- result$Ddnam
    tmp_m <- pmax(beta_val,.Machine$double.eps)/pmax((1-beta_val),.Machine$double.eps)
    m_val <- add_noise_gaussian(log2(tmp_m), sd_dnam)
    Ddnam_noise <- 2^m_val/(2^m_val+1) #probe x sample
    Ddnam_noise[Ddnam_noise<0] <- beta_val[Ddnam_noise<0]
    rownames(Ddnam_noise)=rownames(result$Ddnam)
    colnames(Ddnam_noise)=colnames(result$Ddnam)
    result_noise$Ddnam <- Ddnam_noise
  }
  return(result_noise)
}

add_noise_nb = function(dt, p, sd) {
  delta = matrix(rnorm(prod(dim(dt)), mean=0, sd=sd), nrow=nrow(dt))
  mu_i_0 = t(t(dt)/colSums(dt)) * mean(colSums(dt))
  sigma_i = (1.8*p + 1/sqrt(mu_i_0))*exp(delta/2)
  shape = 1/(sigma_i^2)
  scale = mu_i_0/(shape + .Machine$double.eps)
  mu_i = matrix(rgamma(prod(dim(dt)), shape=shape, scale=scale),
                nrow=nrow(dt))
  v_i = matrix(rpois(prod(dim(dt)), mu_i),
               nrow=nrow(dt))
  return(v_i)
}

add_noise_gaussian = function(dt, sd, mean=0) {
  noise = matrix(rnorm(prod(dim(dt)), mean = mean, sd = sd), nrow = nrow(dt))
  data_noise = dt + noise
  return(data_noise)
}

add_noise_copula <- function(result, empirical, T_rna=NULL, T_dnam=NULL) {
  # compute empirical copula and add noise
  if (!is.null(result$Drna)) {
    N = ncol(result$Drna)
    factor = ncol(result$Drna) / ncol(empirical$Amat)
    epsilon_rna = empirical$Drna - T_rna %*% empirical$Amat
    epsilon_rna_center = do.call(cbind,replicate(factor, epsilon_rna, simplify=FALSE)) - tcrossprod(rowMeans(epsilon_rna), rep(1,N))
    epsilon_rna_scal = epsilon_rna_center/max(1e-8,tcrossprod(sqrt(rowMeans(epsilon_rna_center^2)),rep(1,N)))
    epsilon_rna_pobs = pobs(t(epsilon_rna_scal))
    rownames(epsilon_rna_pobs) = NULL
    fitted_copula_rna = empCopula(epsilon_rna_pobs)
    noise_rna_wo_marg = t(rCopula(N, fitted_copula_rna))

    noise_rna = t(sapply(1:nrow(noise_rna_wo_marg), function(line_i) {
      if (min(epsilon_rna[line_i,])<0) {
        tmp_eps = epsilon_rna[line_i,] - min(epsilon_rna[line_i,])
      } else {tmp_eps = epsilon_rna[line_i,]}
      mean_data = mean(tmp_eps)
      var_data = var(tmp_eps)*(N - 1)/N # unbiased variance
      size_estim = mean_data^2 / max(1e-8,(var_data - mean_data))
      if (size_estim <= 0) { # always with very low data
        size_estim = 0.5
      }
      eps_i = stats::qnbinom(noise_rna_wo_marg[line_i,], size = size_estim, mu = mean_data)
      if (min(epsilon_rna[line_i,])<0) {
        eps_i = eps_i - mean(eps_i) + mean(epsilon_rna[line_i,])
      }
      eps_i}))
    rownames(noise_rna) = rownames(noise_rna_wo_marg)
  
    result$Drna = result$Drna + noise_rna

    which_gene_neg = apply(result$Drna, 1, \(gene_i){any(gene_i<0)})
    which_gene_neg = names(which_gene_neg[which(which_gene_neg)])
    for (gene_i in which_gene_neg) {
      result$Drna[gene_i,] = result$Drna[gene_i,] - min(result$Drna[gene_i,]) + runif(N, min = 0, max = 1e-3)
    }
  }
  if (!is.null(result$Ddnam)) {
    N = ncol(result$Ddnam)
    factor = ncol(result$Ddnam) / ncol(empirical$Amat)
    T_dnam[T_dnam==0]
    bet_val = result$Ddnam
    bet_val_empirical = empirical$Ddnam
    bet_val_empirical[bet_val_empirical==0] = 1e-8
    m_val = log2(bet_val/(1-bet_val))
    m_val_empirical = log2(bet_val_empirical/(1-bet_val_empirical))
    epsilon_met = m_val_empirical - log2(T_dnam/(1-T_dnam)) %*% empirical$Amat
    epsilon_met_center = do.call(cbind,replicate(factor, epsilon_met, simplify=FALSE)) - tcrossprod(rowMeans(epsilon_met), rep(1,N))
    epsilon_met_scal = epsilon_met_center/max(1e-8,tcrossprod(sqrt(rowMeans(epsilon_met_center^2)),rep(1,N)))
    epsilon_met_pobs = pobs(t(epsilon_met_scal))
    rownames(epsilon_met_pobs) = NULL
    fitted_copula_met = empCopula(epsilon_met_pobs)
    noise_met_wo_marg = t(rCopula(N, fitted_copula_met))
    
    noise_met = t(pbapply::pbsapply(1:nrow(noise_met_wo_marg), function(ligne_i) {
      if (min(epsilon_met[ligne_i,])<0) {
        tmp_eps = epsilon_met[ligne_i,] - min(epsilon_met[ligne_i,])
      } else {
        tmp_eps = epsilon_met[ligne_i,]
      }
      mean_data = mean(tmp_eps)
      var_data = var(tmp_eps)*(N - 1)/N # unbiased variance
      eps_i = qnorm(noise_met_wo_marg[ligne_i,], mean = mean_data, sd = var_data)
      if (min(epsilon_met[ligne_i,])<0) {
        eps_i = eps_i - mean(eps_i) + mean(epsilon_met[ligne_i,])
      }
      eps_i}))
    rownames(noise_met) = rownames(noise_met_wo_marg)
    
    m_val = m_val + noise_met
    beta_val <- 2^m_val/(2^m_val+1)
    beta_val[beta_val<0] <- result$Ddnam[beta_val<0]
    beta_val[beta_val>1] <- result$Ddnam[beta_val>1]
    result$Ddnam = beta_val
  }

  return(result)
}

