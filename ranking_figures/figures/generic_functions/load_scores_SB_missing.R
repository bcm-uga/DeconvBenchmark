meth_rna_sup = c("DeconRNASeq", "nnls", "ols","svr","CIBERSORT", "elasticnet", "rlr","WISP", "InstaPrism", "fardeep", "fardeepsto")
meth_rna_unsup = c("ICA", "NMF", "PREDE", "debCAM", "CDSeq")
meth_met_sup = c("rlr","CIBERSORT", "epidishCP","InstaPrism","nnls")
meth_met_unsup = c("RefFreeEWAS", "ICA", "EDec", "MeDeCom", "NMF","debCAM")

## ----
## load scores
## ----
load_data = function(experiment) {
  scores = readRDS(paste0("../1SB_missing/",experiment,"/silico/perf_scores/231027_scores.rds"))
  time = readRDS(paste0("../1SB_missing/",experiment,"/silico/perf_scores/231027_time.rds"))
  
  scores <- scores %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec,"-",experiment),
           score = paste(score,setting),
           dataset = case_when(dataset == "dBREAST" ~ "BrCL1", 
                               dataset == "dPANCREAS" ~ "PaCL1",
                               dataset == "lot1" ~ "PaCL2",
                               dataset == "Cobos" ~ "BrCL2",
                               dataset == "Hoek" ~ "BlCL",
                               dataset == "He" ~ "LuCL"),
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  time <- time %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec,"-",experiment),
           score = paste(score,"time"),
           dataset = case_when(dataset == "dBREAST" ~ "BrCL1", 
                               dataset == "dPANCREAS" ~ "PaCL1",
                               dataset == "lot1" ~ "PaCL2",
                               dataset == "Cobos" ~ "BrCL2",
                               dataset == "Hoek" ~ "BlCL",
                               dataset == "He" ~ "LuCL"),
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  return(list("scores"=scores,"time"=time))
}
