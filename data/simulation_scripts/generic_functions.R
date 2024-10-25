generate_proportions <- function(n_samples, celltypes, alph, varCrit) {
  num.cell.types = length(celltypes)
  # Compute proportions of various compartments
  comp_mix <-t(gtools::rdirichlet(n = n_samples, alpha=alph*varCrit))
  # Merge total matrix proportion
  rownames(comp_mix) = celltypes
  colnames(comp_mix) = paste("train", 1:n_samples, sep = "")
  return(comp_mix)
}

# 0.45 for fib
# 0.1 for imm
# 0.15 for norm
# 0.3 for tum
generate_simu_set <- function(profile_rna, profile_met, n_samples=500, alph=c(0.45, 0.1, 0.15, 0.3), varCrit=10) {
  num.cell.types <- ncol(profile_rna)
  test_solution <- generate_proportions(n_samples, colnames(profile_rna), alph, varCrit) #prop x sample
  test_data_rna <- as.matrix(profile_rna[,colnames(profile_rna)]) %*% test_solution[colnames(profile_rna),] #normalized gene x sample
  test_data_met <- as.matrix(profile_met[,colnames(profile_rna)]) %*% test_solution[colnames(profile_rna),] #probe x sample
  return(list(D_rna=test_data_rna,
              D_met=test_data_met,
              A=test_solution[colnames(profile_rna),]))
}

generate_simu_noise <- function(D_rna, D_met, p, sd_rna=1, sd_met=3) {
  row_names=rownames(D_rna)
  D_rna <- add_noise_nb(D_rna, p, sd_rna)
  beta_val <- D_met
  tmp_m <- pmax(beta_val,.Machine$double.eps)/pmax((1-beta_val),.Machine$double.eps)
  m_val <- add_noise_gaussian(log2(tmp_m), sd_met)
  D_met <- 2^m_val/(2^m_val+1) #probe x sample
  D_met[D_met<0] <- beta_val[D_met<0]
  rownames(D_rna)=row_names
  colnames(D_rna)=colnames(D_met)
  return(list(D_rna=D_rna, D_met=D_met))
}

add_noise_gaussian = function(dt, sd, mean=0) {
  noise = matrix(rnorm(prod(dim(dt)), mean = mean, sd = sd), nrow = nrow(dt))
  data_noise = dt + noise
  return(data_noise)
}

add_noise_nb = function(dt, p, sd=1) {
  delta = matrix(rnorm(prod(dim(dt)), mean = 0, sd = sd),nrow=nrow(dt))
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

scale_data <- function(D_rna, D_met=NULL) {
  size.factor <- DESeq2::estimateSizeFactorsForMatrix(D_rna) 
  D_rna <- sweep(D_rna, 2, size.factor, "/")
  D_rna <- log2(D_rna + 1)
  D_rna <- D_rna/sqrt(nrow(D_rna))
  if (!is.null(D_met)) {
    D_met <- t(scale(t(D_met)))
    D_met <- D_met/sqrt(nrow(D_met))
    return(list(D_rna=D_rna, D_met=D_met))
  }
  else return(D_rna)
}
