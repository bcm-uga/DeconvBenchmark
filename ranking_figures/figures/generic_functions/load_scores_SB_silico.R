## ----
## load scores
## ----
load_data = function(date, score_path) {
  scores = readRDS(paste0(score_path,date,"_silico_scores.rds"))
  time = readRDS(paste0(score_path,date,"_silico_time.rds"))
  
  scores <- scores %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec),
           score = paste(score,setting),
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  time <- time %>%
    mutate(candidate = paste0(block,"-",deconv,"-",feat_selec),
           score = paste(score,"time"),
           dataset = paste0(dataset,'-',block)) %>%
    select(values,score,sim,dataset,candidate)
  return(list("scores"=scores,"time"=time))
}
