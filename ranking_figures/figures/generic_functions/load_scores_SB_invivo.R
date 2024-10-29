meth_rna_sup = c("DeconRNASeq", "nnls", "ols","svr","CIBERSORT", "elasticnet", "rlr","WISP", "InstaPrism", "fardeep", "fardeepsto")
meth_rna_unsup = c("ICA", "NMF", "PREDE", "debCAM", "CDSeq")
meth_met_sup = c("rlr","CIBERSORT", "epidishCP","InstaPrism","nnls")
meth_met_unsup = c("RefFreeEWAS", "ICA", "EDec", "MeDeCom", "NMF","debCAM")

## ----
## load scores
## ----
load_data = function() {
  scores = readRDS("../1SB_invivo/perf_scores/scores.rds")
  time = readRDS("../1SB_invivo/perf_scores/time.rds")
  
  scores <- scores %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec),
           score = paste(score,setting),
           sim = 'na',
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  time <- time %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec),
           score = paste(score,"time"),
           sim = 'na',
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  return(list("scores"=scores,"time"=time))
}
