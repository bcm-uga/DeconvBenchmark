library(dplyr)
library(matrixStats)


function_metrics <- function(real, prediction, method = c("rmse","mae","pearson")) {
  method <- match.arg(method)
  prediction[is.na(prediction)] = 0
  if (method=="rmse") {
    return(sqrt(mean((c(as.matrix(real)) - c(as.matrix(prediction)))^2)))
  }
  else if (method=="mae") {
    return(mean(abs(c(as.matrix(real)) - c(as.matrix(prediction)))))
  }
  else if (method=="pearson") {
    return(cor(c(as.matrix(real)), c(as.matrix(prediction)), method="pearson"))
  }
}

compute_scores <- function (A_true, A_pred, score_methods=c("rmse", "mae", "pearson")) {
  if (!all(is.na(A_pred))) {
    if (nrow(A_pred) != nrow(A_true)) {
      warning("The number of predicted cell types doesn't match the real number of cell types")
      if (nrow(A_pred)==nrow(A_true)-1) {
        print("There is one less cell type in A_pred")
        A_true <- A_true[sort(rownames(A_pred)),]
        A_pred <- A_pred[sort(rownames(A_pred)),]
      }
      else if (nrow(A_pred)==nrow(A_true)+1) {
        print("There is one more cell type in A_pred")
        A_true <- A_true[sort(rownames(A_true)),]
        A_pred <- A_pred[sort(rownames(A_true)),]
      }
      else {stop("Error in the number of predicted cell types")}
    } else {A_pred <- A_pred[rownames(A_true), ]}

    score1 <- sapply(score_methods, function(x) function_metrics(A_true, A_pred, x)) #perf_g
    score_celltype <- sapply(seq(nrow(A_pred)), function(i) function_metrics(A_true[i,], A_pred[i,], "pearson"))
    score_celltype_sd <- sd(score_celltype, na.rm=T) #sd_c
    score_celltype_median <- median(score_celltype, na.rm=T) #med_c
    score_sample <- sapply(seq(ncol(A_pred)), function(i) function_metrics(A_true[, i], A_pred[, i], "pearson"))
    score_sample_sd <- sd(score_sample, na.rm=T) #sd_s
    score_sample_median <- median(score_sample, na.rm=T) #med_s
  } else {
    score1 <- rep(NA, length(score_methods))
    names(score1) <- score_methods
    score_celltype_sd <- NA
    score_celltype_median <- NA
    score_sample_sd <- NA
    score_sample_median <- NA
  }

  return(data.frame("values" = c(score1,
                                 score_celltype_sd, score_celltype_median,
                                 score_sample_sd, score_sample_median),
                    "score" = c(score_methods,
                                rep("pearson",4)),
                    "setting" = c(rep("perf_g",length(score_methods)), "sd_c", "med_c", "sd_s", "med_s")))
}

setting_SB_silico <- function(data_path, deconv_path, score_path, date, score_methods=c("rmse", "mae", "pearson"), fs_bool=TRUE) {
  blocks <- list.dirs(deconv_path, recursive = F, full.names = F)
  meth_classes <- list.dirs(paste0(deconv_path,blocks[1]), recursive = F, full.names = F)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    print(paste0("block ",block))
    datasets <- unique(sapply(list.files(paste0(data_path,block)), function(x)
      strsplit(x, "_")[[1]][2]))
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      print(paste0("class ",meth_class))
      # for all datasets
      df_class[[meth_class]] <- do.call(rbind,lapply(datasets, function(Data) {
        print(paste0("**  ",Data,"  **"))
        sim_files <- list.files(paste0(data_path,block), pattern = paste0(date, "_", Data))
        # for all simulations
        df_lot <- do.call(rbind,lapply(sim_files, function(Sim_file) {
          sim <- strsplit(strsplit(Sim_file, ".rds")[[1]], "_sim")[[1]][[2]]
          A_true <- readRDS(paste0(data_path,block,"/", Sim_file))$A_ref
          # for all deconvolution methods
          res_files <- list.files(paste0(deconv_path,block,'/',meth_class, "/"),
                                  pattern = glob2rx(paste0(date, "_", Data, "_Apred_", "*_*", "_sim", sim, ".rds")))
          methods = unique(sapply(res_files, function(x) strsplit(x,"_")[[1]][5]))
          df_sim <- do.call(rbind,lapply(methods, function (Method) {
            fs = unique(sapply(res_files, function(x) strsplit(x,"_")[[1]][4]))
            df_meth <- do.call(rbind,lapply(fs, function(i) {
              res_file = res_files[intersect(grep(i,res_files),
                                             grep(Method,res_files))]
              A_pred <- readRDS(paste0(deconv_path,block,'/', meth_class, "/",
                                       date, "_", Data, "_Apred_", i, "_", Method, "_sim", sim, ".rds"))
              if (Data=='PaCL2') {
                rownames(A_pred)=gsub("\\."," ",rownames(A_pred))
                rownames(A_pred)=gsub("_"," ",rownames(A_pred))
                }
              df_fs <- compute_scores(A_true, A_pred, score_methods)
              df_fs$feat_selec <- i
              return(df_fs)
            }))
            df_meth$deconv <- Method
            return(df_meth)
          }))
          df_sim$sim <- as.integer(sim)
          return(df_sim)
        }))
        df_lot$dataset <- Data
        return(df_lot)
      }))
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  score_perf_df = do.call(rbind,df_block)
  saveRDS(score_perf_df, paste0(score_path,date, "_scores.rds"))
}

generate_score_SB_mean <- function(data_path, deconv_path, date, score_methods, deconv_methods) {
  path_suffix = deconv_methods %>% select(Block,Supervised) %>% unique()
  lots <- lapply(deconv_methods %>%
                   filter(!duplicated(paste(Block,Supervised))) %>%
                   pull(Block), function(block)
    unique(sapply(list.files(paste0(data_path, block,"/")), function(file)
      strsplit(file,"_")[[1]][2])))
  
  # for all path_suffixes
  score_perf_setting <- do.call(rbind,lapply(seq(nrow(path_suffix)), function (path) {
    block = path_suffix$Block[path]
    class = path_suffix$Supervised[path]
    print(paste("Running setting",block,class))
    # for all lots
    score_perf_lot = lapply(lots[[path]], function(lot) {
      print(paste("Running lot",lot))
      sim_files <- list.files(paste0(data_path,block), pattern = paste0(date, "_", lot,"_sim"))
      # for all simulations
      score_perf_sim <- lapply(sim_files, function(sim_file) {
        sim <- strsplit(strsplit(sim_file, ".rds")[[1]], "_sim")[[1]][[2]]
        A_true <- readRDS(paste0(paste0(data_path,block,"/"), sim_file))$A_ref
        # in 1 class
        list_meth = deconv_methods %>%
          filter(Block==block,
                 Supervised==class) %>%
          select(DeconvTool,FS)
        A_pred <- list()
        for (i in seq(nrow(list_meth))) {
          A_pred[[list_meth$DeconvTool[i]]] <- readRDS(paste0(deconv_path, block, "/",class,"/", date,"/", date, "_", lot, "_Apred",list_meth$FS[i],"_", list_meth$DeconvTool[i], "_sim", sim, ".rds"))
        }
        A_pred[[length(A_pred)+1]] = Reduce("+",A_pred)/length(A_pred)
        names(A_pred) <- c(list_meth$DeconvTool,paste("Consensus",class))
        df_inter <- lapply(A_pred, function(x)
          compute_scores(A_true, x, score_methods))
        df_scores = do.call(rbind,mapply(function(x,y) {
          x$candidate = y
          x
        }, df_inter, names(A_pred), SIMPLIFY=F))
        df_scores$sim <- as.integer(sim)
        return(df_scores)
      })
      score_perf_sim <- do.call(rbind, score_perf_sim)
      score_perf_sim$dataset <- lot
      return(score_perf_sim)
      })
    score_perf_lot <- do.call(rbind, score_perf_lot) %>%
      mutate(score = paste(score, setting),
             block = block,
             class = class) %>%
      select(values, score, sim, dataset, candidate, block, class)
    score_perf_lot
  }))
  return(score_perf_setting)
}

generate_score_SB_invitro <- function(data_path, deconv_path, perf_score_path, score_methods) {
  blocks <- list.dirs(deconv_path, recursive = F, full.names = F)
  meth_classes <- list.dirs(paste0(deconv_path,blocks[1]), recursive = F, full.names = F)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    print(paste0("block ",block))
    lots <- list.files(data_path)
    lots = lots[grep(block,lots)]
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      print(paste0("class ",meth_class))
      # for all datasets
      df_class[[meth_class]] <- do.call(rbind,lapply(lots, function(Data) {
        A_true <- readRDS(paste0(data_path,Data,"/Amat.rds"))
        # for all deconvolution methods
        res_files <- list.files(paste0(deconv_path,block,'/',meth_class, "/"),
                                  pattern = strsplit(Data,"_")[[1]][1])
        methods = unique(sapply(res_files, function(x) strsplit(strsplit(x,"_")[[1]][3],".rds")[[1]][1]))
        df_lot <- do.call(rbind,lapply(methods, function (Method) {
          # for all feature selection strategies, if applicable
          print(Method)
          fs = c("none","toast","hvf")
          df_meth <- do.call(rbind,lapply(fs, function(i) {
            res_file = res_files[intersect(grep(i,res_files),
                                             grep(paste0(Method,".rds"),res_files))]
            A_pred <- readRDS(paste0(deconv_path,block,'/', meth_class, "/", res_file))
            df_fs <- compute_scores(A_true, A_pred, score_methods)
            df_fs$feat_selec <- i
            df_fs
            }))
          df_meth$deconv <- Method
          return(df_meth)
        }))
        df_lot$dataset <- Data
        return(df_lot)
      }))
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  score_perf_df = do.call(rbind,df_block)
  saveRDS(score_perf_df, paste0(perf_score_path, "scores.rds"))
}

generate_score_SB_invivo <- function(data_path, deconv_path, perf_score_path, score_methods) {
  blocks <- list.dirs(deconv_path, recursive = F, full.names = F)
  meth_classes <- list.dirs(paste0(deconv_path,blocks[1]), recursive = F, full.names = F)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    print(paste0("block ",block))
    lots <- list.files(data_path)
    lots = lots[grep(block,lots)]
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      print(paste0("class ",meth_class))
      # for all datasets
      df_class[[meth_class]] <- do.call(rbind,lapply(lots, function(Data) {
        A_true <- readRDS(paste0(data_path,Data,"/Amat.rds"))
        # for all deconvolution methods
        res_files <- list.files(paste0(deconv_path,block,'/',meth_class, "/"),
                                pattern = strsplit(Data,"_")[[1]][1])
        methods = unique(sapply(res_files, function(x) strsplit(strsplit(x,"_")[[1]][3],".rds")[[1]][1]))
        df_lot <- do.call(rbind,lapply(methods, function (Method) {
          # for all feature selection strategies, if applicable
          print(Method)
          fs = c("none","toast","hvf")
          df_meth <- do.call(rbind,lapply(fs, function(i) {
            res_file = res_files[intersect(grep(i,res_files),
                                           grep(paste0(Method,".rds"),res_files))]
            A_pred <- readRDS(paste0(deconv_path,block,'/', meth_class, "/", res_file))
            if (Data=="SkREAL_rna" & meth_class=="sup") {
              OtherCells = colSums(A_pred[grep("Other",rownames(A_pred)),])
              A_pred = rbind(A_pred[-grep("Other",rownames(A_pred)),],
                             OtherCells)
            }
            df_fs <- compute_scores(A_true, A_pred, score_methods)
            df_fs$feat_selec <- i
            df_fs
          }))
          df_meth$deconv <- Method
          return(df_meth)
        }))
        df_lot$dataset <- Data
        return(df_lot)
      }))
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  score_perf_df = do.call(rbind,df_block)
  saveRDS(score_perf_df, paste0(perf_score_path, "scores.rds"))
}

generate_time_SB_silico <- function(timing_path, perf_score_path, date) {
  blocks <- list.dirs(timing_path, recursive = FALSE, full.names = FALSE)
  meth_classes <- list.dirs(paste0(timing_path,blocks[1]), recursive = FALSE, full.names = FALSE)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      df_class[[meth_class]] <- do.call(rbind,lapply(list.files(paste0(timing_path,block,'/', meth_class,"/",date),
                                                                pattern = glob2rx(paste0("*_timing*.rds"))),
                                                     function (time_file) {
                                                       suffix <- tail(strsplit(strsplit(time_file, "_")[[1]][3], "timing")[[1]], 1)
                                                       lot <- strsplit(time_file, "_")[[1]][2]
                                                       sim <- as.numeric(strsplit(strsplit(time_file, "_sim")[[1]][2], ".rds")[[1]][1])
                                                       meth <- strsplit(time_file, "_")[[1]][4]
                                                       time <- readRDS(paste0(timing_path,block,'/', meth_class, "/", date, "/", time_file))
                                                       if (length(grep("LessSamples",date))==0) {
                                                         time = time/120
                                                       } else {time = time/30}
                                                       if (is.null(time)) time <- NA
                                                       res_df <- data.frame("values" = time,
                                                                            "score" = "time",
                                                                            "deconv" = meth,
                                                                            "dataset" = lot,
                                                                            "sim" = sim,
                                                                            "feat_selec" = suffix)
                                                       return(res_df)
                                                     }))
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  file_path <- paste0(perf_score_path, date, "_time.rds")
  saveRDS(do.call(rbind,df_block), file_path)
}

generate_time_SB_mean <- function(data_path, timing_path, date, deconv_methods) {
  path_suffix = deconv_methods %>% select(Block,Supervised) %>% unique()
  lots <- lapply(deconv_methods %>%
                   filter(!duplicated(paste(Block,Supervised))) %>%
                   pull(Block), function(block)
                     unique(sapply(list.files(paste0(data_path, block,"/")), function(file)
                       strsplit(file,"_")[[1]][2])))
  
  # for all path_suffixes
  time_setting <- do.call(rbind,lapply(seq(nrow(path_suffix)), function (path) {
    block = path_suffix$Block[path]
    class = path_suffix$Supervised[path]
    print(paste("Running setting",block,class))
    # for all lots
    time_lot = do.call(rbind,lapply(lots[[path]], function(lot) {
      sim_files <- list.files(paste0(data_path,block), pattern = paste0(date, "_", lot,"_sim"))
      # for all simulations
      time_sim <- lapply(sim_files, function(sim_file) {
        sim <- strsplit(strsplit(sim_file, ".rds")[[1]], "_sim")[[1]][[2]]
        # in 1 class, all methods
        list_meth = deconv_methods %>%
          filter(Block==block,
                 Supervised==class) %>%
          select(DeconvTool,FS)
        time_res = rep(NA,nrow(list_meth)+1)
        for (i in seq(nrow(list_meth))) {
          time_res[i] <- readRDS(paste0(timing_path, block, "/",class,"/", date,"/", date, "_", lot, "_timing",list_meth$FS[i],"_", list_meth$DeconvTool[i], "_sim", sim, ".rds"))
        }
        time_res[length(time_res)] = sum(time_res[seq(length(time_res)-1)])
        if (length(grep("LessSamples",date))==0) {
          time_res = time_res/120
        } else {time_res = time_res/30}
        data.frame(values=time_res,
                   score='time time',
                   sim=as.numeric(sim),
                   dataset=lot,
                   candidate=c(list_meth$DeconvTool,paste0("Consensus ",class)))
      })
      return(do.call(rbind,time_sim))
    }))
    time_lot$block=block
    time_lot$class=class
    time_lot
  }))
  return(time_setting)
}

generate_time_invitro <- function(timing_path, perf_score_path, input_data_path="data/") {
  blocks <- list.dirs(timing_path, recursive = FALSE, full.names = FALSE)
  meth_classes <- list.dirs(paste0(timing_path,blocks[1]), recursive = FALSE, full.names = FALSE)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      df_class[[meth_class]] <- do.call(rbind,lapply(list.files(paste0(timing_path, block, "/", meth_class),
                                                                pattern = glob2rx(paste0("*_timing*.rds"))),
                                                     function (time_file) {
                                                       suffix <- strsplit(strsplit(time_file, "timing")[[1]][2],"_")[[1]][1]
                                                       lot <- strsplit(time_file, "_")[[1]][1]
                                                       lot <- grep(lot,list.files(input_data_path), value=T)
                                                       meth <- strsplit(strsplit(time_file, "_")[[1]][3], ".rds")[[1]][1]
                                                       time <- readRDS(paste0(timing_path, block, "/", meth_class, "/", time_file))
                                                       n_sample = ifelse(lot=="BlMIX_met",12,
                                                                         ifelse(lot=="PaMIX_met_rna",30,
                                                                                ifelse(lot=="BrMIX_rna",18,warning("issue with lot"))))
                                                       time = time/n_sample
                                                       if (is.null(time)) time <- NA
                                                       return(data.frame(values = time,
                                                                         score = "time",
                                                                         deconv = meth,
                                                                         dataset = lot,
                                                                         feat_selec = suffix,
                                                                         class=meth_class))
                                                     }))
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  file_path <- paste0(perf_score_path, "time.rds")
  saveRDS(do.call(rbind,df_block), file_path)
}

generate_time_invivo <- function(timing_path, perf_score_path) {
  blocks <- list.dirs(timing_path, recursive = FALSE, full.names = FALSE)
  meth_classes <- list.dirs(paste0(timing_path,blocks[1]), recursive = FALSE, full.names = FALSE)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      df_class[[meth_class]] <- do.call(rbind,lapply(list.files(paste0(timing_path, block, "/", meth_class),
                                                                pattern = glob2rx(paste0("*_timing*.rds"))),
                                                     function (time_file) {
                                                       suffix <- strsplit(strsplit(time_file, "timing")[[1]][2],"_")[[1]][1]
                                                       lot <- strsplit(time_file, "_")[[1]][1]
                                                       lot <- grep(lot,list.files("data/"), value=T)
                                                       meth <- strsplit(strsplit(time_file, "_")[[1]][3], ".rds")[[1]][1]
                                                       time <- readRDS(paste0(timing_path, block, "/", meth_class, "/", time_file))
                                                       n_sample = ifelse(lot=="BlREAL1_rna",5,
                                                                         ifelse(lot=="BlREAL2_met",12,
                                                                                ifelse(lot=="BlREAL3_met",335,
                                                                                       ifelse(lot=="SkREAL_rna",4,warning("issue with lot")))))
                                                       time = time/n_sample
                                                       if (is.null(time)) time <- NA
                                                       return(data.frame(values = time,
                                                                         score = "time",
                                                                         deconv = meth,
                                                                         dataset = lot,
                                                                         feat_selec = suffix,
                                                                         class=meth_class))
                                                     }))
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  file_path <- paste0(perf_score_path, "time.rds")
  saveRDS(do.call(rbind,df_block), file_path)
}