## ----
## load scores
## ----
load_data_missing = function(experiment, score_path) {
  scores = readRDS(paste0(score_path,experiment,"_scores.rds"))
  time = readRDS(paste0(score_path,experiment,"_time.rds"))
  
  scores <- scores %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec,"-",experiment),
           score = paste(score,setting),
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  time <- time %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec,"-",experiment),
           score = paste(score,"time"),
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  return(list("scores"=scores,"time"=time))
}
