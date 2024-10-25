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
  mu_i_0 = dt/colSums(dt) * mean(colSums(dt))
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
