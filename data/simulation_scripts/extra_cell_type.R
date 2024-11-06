set.seed(0)

## ----
## Parameters, put your own
## ----
date="241025"
rna_data = c("BlCL","BrCL1","BrCL2","PaCL1","PaCL2")
  
## ----
## load original Tref
## ----
Tref = lapply(list.files("../references", full.names = T, pattern = "CL"),
              readRDS)
names(Tref) = list.files("../references", pattern = "CL")


## ----
## Functions
## ----
add_noise_gaussian = function(dt, sd, mean=0) {
  noise = matrix(rnorm(prod(dim(dt)), mean = mean, sd = sd), nrow = nrow(dt))
  data_noise = dt + noise
  return(data_noise)
}

add_noise_nb = function(dt, sd, p=.1) {
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

add_noise <- function(mat, block, sd_rna=1, sd_met=1) {
  if (block=='rna') {
    mat_noisy <- add_noise_nb(mat, sd_rna)
  }
  if (block=='dnam') {
    m_val <- add_noise_gaussian(log2(mat/(1-mat)), sd_met)
    mat_noisy <- 2^m_val/(2^m_val+1) #probe x sample
    mat_noisy[mat_noisy<0] <- mat[mat_noisy<0]
    
  }
  rownames(mat_noisy) = rownames(mat)
  colnames(mat_noisy) = colnames(mat)
  return(mat_noisy)
}

## ----
## add fake cell types generated based on immune cell types, adapt if new datasets
## ----
sapply(Tref,colnames)
imm_cell_types = list(c('B','mDC','Mono','Neut','NK','T'),
                      'lymphocytes',
                      'lymphocytes',
                      c('Jurkat','Thp1'),
                      c('CD19B','CD4T','CD8T','Monocyte','Neutrophil','NKcell','WBC'),
                      'IMMUNE',
                      'IMMUNE',
                      c('B.cells','CD4.T.cells','CD8.T.cells','Macrophages','Neutrophils'),
                      c('B.cells','CD4.T.cells','CD8.T.cells','Macrophages','Neutrophils'))

for (data in seq_along(Tref)) {
  name_data = gsub(".rds","",names(Tref)[data])
  print(name_data)
  block = ifelse(name_data %in% rna_data,"rna","dnam")
  Tmat_imm = Tref[[data]][,imm_cell_types[[data]]]
  if (length(imm_cell_types[[data]])==1) {
    Tmat_imm_meta = add_noise(matrix(Tmat_imm),block=block)
  } else {
    Tmat_imm_meta = add_noise(matrix(rowMeans(Tmat_imm)),block=block)
  }
  Tmat_added = cbind(Tref[[data]],Tmat_imm_meta)
  colnames(Tmat_added) = c(colnames(Tref[[data]]),"MetaImmune")
  saveRDS(Tmat_added,
          paste0("../references/", name_data, "_added.rds"))
  
}
