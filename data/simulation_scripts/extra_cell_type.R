set.seed(0)

## ----
## List all silico datasets
## ----
blocks = list.files("../../../acacia_2final/results/0simu/simulations", full.names = T)
lots = lapply(blocks, function(x) unique(sapply(list.files(x, pattern="231027"), function(y)
  strsplit(y,'_')[[1]][2])))
blocks = list.files("../../../acacia_2final/results/0simu/simulations")
names(lots) = blocks

## ----
## load original Tref
## ----
Tref = lapply(names(lots), function(block)
  lapply(lots[[block]], function(lot)
    data.frame(readRDS(paste0("../../../acacia_2final/results/0simu/simulations/",block,"/231027_", lot, "_T_",block,"_ref.rds")))))
names(Tref) = blocks

## ----
## add fake cell types generated based on immune cell types
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
  if (block=='met') {
    m_val <- add_noise_gaussian(log2(mat/(1-mat)), sd_met)
    mat_noisy <- 2^m_val/(2^m_val+1) #probe x sample
    mat_noisy[mat_noisy<0] <- mat[mat_noisy<0]
    
  }
  rownames(mat_noisy) = rownames(mat)
  colnames(mat_noisy) = colnames(mat)
  return(mat_noisy)
}

lapply(Tref, function(x) sapply(x,colnames))
imm_cell_types   = list('met'=list('lymphocytes',
                                   'IMMUNE',
                                   c('CD19B','CD4T','CD8T','Monocyte','Neutrophil','NKcell','WBC'),
                                   c('B.cells','CD4.T.cells','CD8.T.cells','Macrophages','Neutrophils')),
                        'rna'=list(c('Jurkat','Thp1'),
                                   'lymphocytes',
                                   'IMMUNE',
                                   c('B','mDC','Mono','Neut','NK','T'),
                                   c('B.cells','CD4.T.cells','CD8.T.cells','Macrophages','Neutrophils')))

for (block in blocks) {
  for (dataset in seq_along(Tref[[block]])) {
    print(lots[[block]][[dataset]])
    Tmat_imm = Tref[[block]][[dataset]][,imm_cell_types[[block]][[dataset]]]
    if (length(imm_cell_types[[block]][[dataset]])==1) {
      Tmat_imm_meta = add_noise(matrix(Tmat_imm),block=block)
    } else {
      Tmat_imm_meta = add_noise(matrix(rowMeans(Tmat_imm)),block=block)
    }
    Tmat_addedinT = cbind(Tref[[block]][[dataset]],Tmat_imm_meta)
    colnames(Tmat_addedinT) = c(colnames(Tref[[block]][[dataset]]),"MetaImmune")
    saveRDS(Tmat_addedinT,
            paste0("../../../acacia_2final/results/0simu/simulations/",block,"/231027_", lots[[block]][[dataset]], "_T_",block,"_ref_addedinT.rds"))
  }
}
