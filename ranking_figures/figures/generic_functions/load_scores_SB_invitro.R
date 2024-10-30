## ----
## load scores
## ----
load_data_vitro = function(score_path) {
  scores = readRDS(paste0(score_path,"vitro_scores.rds"))
  time = readRDS(paste0(score_path,"vitro_time.rds"))
  
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
