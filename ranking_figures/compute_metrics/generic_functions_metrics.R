library(dplyr)
library(matrixStats)


metrics <- function(real, prediction, method = c("rmse","mae","pearson")) {
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

    score1 <- sapply(score_methods, function(x) metrics(A_true, A_pred, x)) #perf_g
    score_celltype <- sapply(seq(nrow(A_pred)), function(i) metrics(A_true[i,], A_pred[i,], "pearson"))
    score_celltype_sd <- sd(score_celltype, na.rm=T) #sd_c
    score_celltype_median <- median(score_celltype, na.rm=T) #med_c
    score_sample <- sapply(seq(ncol(A_pred)), function(i) metrics(A_true[, i], A_pred[, i], "pearson"))
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

score_SB_silico <- function(data_path, deconv_path, score_path, date, score_methods=c("rmse", "mae", "pearson"), fs_bool=TRUE) {
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
      df_class[[meth_class]] <- lapply(datasets, function(Data) {
        print(paste0("**  ",Data,"  **"))
        sim_files <- list.files(paste0(data_path,block), pattern = paste0(date, "_", Data))
        # for all simulations
        df_dataset <- do.call(rbind,lapply(sim_files, function(Sim_file) {
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
        df_dataset$dataset <- Data
        return(df_dataset)
      })
      df_class = df_class[!is.null(df_class)]
      df_class = do.call(rbind,df_class)
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  score_perf_df = do.call(rbind,df_block)
  saveRDS(score_perf_df, paste0(score_path,date, "_silico_scores.rds"))
}

score_SB_consensus <- function(data_path, deconv_path, score_path, date, deconv_methods, score_methods) {
  path_suffix = deconv_methods %>% select(Block,Supervised) %>% unique()
  datasets <- lapply(deconv_methods %>%
                   filter(!duplicated(paste(Block,Supervised))) %>%
                   pull(Block), function(block)
    unique(sapply(list.files(paste0(data_path, block,"/")), function(file)
      strsplit(file,"_")[[1]][2])))
  
  # for all path_suffixes
  score_perf_setting <- do.call(rbind,lapply(seq(nrow(path_suffix)), function (path) {
    block = path_suffix$Block[path]
    class = path_suffix$Supervised[path]
    print(paste("Running setting",block,class))
    # for all datasets
    score_perf_data = lapply(datasets[[path]], function(Data) {
      print(paste("Running dataset",Data))
      sim_files <- list.files(paste0(data_path,block), pattern = paste0(date, "_", Data,"_sim"))
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
          A_pred[[list_meth$DeconvTool[i]]] <- readRDS(paste0(deconv_path, block, "/",class,"/", date, "_", Data, "_Apred_",list_meth$FS[i],"_", list_meth$DeconvTool[i], "_sim", sim, ".rds"))
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
      score_perf_sim$dataset <- Data
      return(score_perf_sim)
      })
    score_perf_data <- do.call(rbind, score_perf_data) %>%
      mutate(score = paste(score, setting),
             block = block,
             class = class) %>%
      select(values, score, sim, dataset, candidate, block, class)
    score_perf_data
  }))
  saveRDS(score_perf_setting, paste0(score_path,date, "_silico_scores_consensus.rds"))
}

score_SB_invitro <- function(data_path, deconv_path, score_path, score_methods) {
  blocks <- list.dirs(deconv_path, recursive = F, full.names = F)
  meth_classes <- list.dirs(paste0(deconv_path,blocks[1]), recursive = F, full.names = F)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    print(paste0("block ",block))
    datasets <- list.files(data_path, pattern="_D_")
    datasets = sapply(datasets[grep(block,datasets)], function(x) strsplit(x,"_")[[1]][1])
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      print(paste0("class ",meth_class))
      # for all datasets
      df_class[[meth_class]] <- do.call(rbind,lapply(datasets, function(Data) {
        A_true <- readRDS(paste0(data_path,Data,"_A.rds"))
        # for all deconvolution methods
        res_files <- list.files(paste0(deconv_path,block,'/',meth_class, "/"),
                                  pattern = Data)
        methods = unique(sapply(res_files, function(x) strsplit(strsplit(x,"_")[[1]][4],".rds")[[1]][1]))
        df_data <- do.call(rbind,lapply(methods, function (Method) {
          # for all feature selection strategies, if applicable
          print(Method)
          fs = unique(sapply(res_files, function(x) strsplit(x,"_")[[1]][3]))
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
        df_data$dataset <- Data
        return(df_data)
      }))
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  score_perf_df = do.call(rbind,df_block)
  saveRDS(score_perf_df, paste0(score_path, "vitro_scores.rds"))
}

score_SB_invivo <- function(data_path, deconv_path, score_path, score_methods) {
  blocks <- list.dirs(deconv_path, recursive = F, full.names = F)
  meth_classes <- list.dirs(paste0(deconv_path,blocks[1]), recursive = F, full.names = F)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    print(paste0("block ",block))
    datasets <- list.files(data_path, pattern="_D_")
    datasets = sapply(datasets[grep(block,datasets)], function(x) strsplit(x,"_")[[1]][1])
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      print(paste0("class ",meth_class))
      # for all datasets
      df_class[[meth_class]] <- do.call(rbind,lapply(datasets, function(Data) {
        A_true <- readRDS(paste0(data_path,Data,"_A.rds"))
        # for all deconvolution methods
        res_files <- list.files(paste0(deconv_path,block,'/',meth_class, "/"),
                                pattern = Data)
        methods = unique(sapply(res_files, function(x) strsplit(strsplit(x,"_")[[1]][4],".rds")[[1]][1]))
        df_dataset <- do.call(rbind,lapply(methods, function (Method) {
          # for all feature selection strategies, if applicable
          print(Method)
          fs = unique(sapply(res_files, function(x) strsplit(x,"_")[[1]][3]))
          df_meth <- do.call(rbind,lapply(fs, function(i) {
            res_file = res_files[intersect(grep(i,res_files),
                                           grep(paste0(Method,".rds"),res_files))]
            A_pred <- readRDS(paste0(deconv_path,block,'/', meth_class, "/", res_file))
            if (Data=="SkREAL" & meth_class=="sup") {
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
        df_dataset$dataset <- Data
        return(df_dataset)
      }))
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  score_perf_df = do.call(rbind,df_block)
  saveRDS(score_perf_df, paste0(score_path, "vivo_scores.rds"))
}

time_SB_silico <- function(timing_path, score_path, date, n_sample=120) {
  blocks <- list.dirs(timing_path, recursive = FALSE, full.names = FALSE)
  meth_classes <- list.dirs(paste0(timing_path,blocks[1]), recursive = FALSE, full.names = FALSE)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      df_class[[meth_class]] <- do.call(rbind,lapply(list.files(paste0(timing_path,block,'/', meth_class),
                                                                pattern = "CL"),
                                                     function (time_file) {
                                                       fs <- strsplit(time_file, "_")[[1]][4]
                                                       data <- strsplit(time_file, "_")[[1]][2]
                                                       sim <- as.numeric(strsplit(strsplit(time_file, "_sim")[[1]][2], ".rds")[[1]][1])
                                                       method <- strsplit(time_file, "_")[[1]][5]
                                                       time <- readRDS(paste0(timing_path,block,'/', meth_class, "/", time_file))
                                                       time = time/n_sample
                                                       if (is.null(time)) time <- NA
                                                       res_df <- data.frame("values" = time,
                                                                            "score" = "time",
                                                                            "deconv" = method,
                                                                            "dataset" = data,
                                                                            "sim" = sim,
                                                                            "feat_selec" = fs)
                                                       return(res_df)
                                                     }))
      df_class[[meth_class]]$class = meth_class
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  saveRDS(do.call(rbind,df_block), paste0(score_path, date, "_silico_time.rds"))
}

time_SB_consensus <- function(data_path, timing_path, score_path, date, deconv_methods, n_sample=120) {
  path_suffix = deconv_methods %>% select(Block,Supervised) %>% unique()
  datasets <- lapply(deconv_methods %>%
                   filter(!duplicated(paste(Block,Supervised))) %>%
                   pull(Block), function(block)
                     unique(sapply(list.files(paste0(data_path, block,"/")), function(file)
                       strsplit(file,"_")[[1]][2])))
  # for all path_suffixes
  time_setting <- do.call(rbind,lapply(seq(nrow(path_suffix)), function (path) {
    block = path_suffix$Block[path]
    class = path_suffix$Supervised[path]
    print(paste("Running setting",block,class))
    # for all datasets
    time_dataset = do.call(rbind,lapply(datasets[[path]], function(Data) {
      sim_files <- list.files(paste0(data_path,block), pattern = paste0(date, "_", Data,"_sim"))
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
          time_res[i] <- readRDS(paste0(timing_path, block, "/",class,"/", date, "_", Data, "_timing_",list_meth$FS[i],"_", list_meth$DeconvTool[i], "_sim", sim, ".rds"))
        }
        time_res[length(time_res)] = sum(time_res[seq(length(time_res)-1)])
        time_res = time_res/n_sample
        data.frame(values=time_res,
                   score='time time',
                   sim=as.numeric(sim),
                   dataset=Data,
                   candidate=c(list_meth$DeconvTool,paste0("Consensus ",class)))
      })
      return(do.call(rbind,time_sim))
    }))
    time_dataset$block=block
    time_dataset$class=class
    time_dataset
  }))
  saveRDS(time_setting, paste0(score_path, date, "_silico_time_consensus.rds"))
}

time_SB_invitro <- function(timing_path, score_path) {
  blocks <- list.dirs(timing_path, recursive = FALSE, full.names = FALSE)
  meth_classes <- list.dirs(paste0(timing_path,blocks[1]), recursive = FALSE, full.names = FALSE)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      df_class[[meth_class]] <- do.call(rbind,lapply(list.files(paste0(timing_path, block, "/", meth_class),
                                                                pattern='MIX'),
                                                     function (time_file) {
                                                       fs <- strsplit(time_file, "_")[[1]][3]
                                                       data <- strsplit(time_file, "_")[[1]][1]
                                                       meth <- strsplit(strsplit(time_file, "_")[[1]][4], ".rds")[[1]][1]
                                                       time <- readRDS(paste0(timing_path, block, "/", meth_class, "/", time_file))
                                                       n_sample = ifelse(data=="BlMIX",12,
                                                                         ifelse(data=="PaMIX",30,
                                                                                ifelse(data=="BrMIX",18,warning("You should specify the number of samples of your dataset."))))
                                                       time = time/n_sample
                                                       if (is.null(time)) time <- NA
                                                       return(data.frame(values = time,
                                                                         score = "time",
                                                                         deconv = meth,
                                                                         dataset = data,
                                                                         feat_selec = fs,
                                                                         class=meth_class))
                                                     }))
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  saveRDS(do.call(rbind,df_block), paste0(score_path, "vitro_time.rds"))
}

time_SB_invivo <- function(timing_path, score_path) {
  blocks <- list.dirs(timing_path, recursive = FALSE, full.names = FALSE)
  meth_classes <- list.dirs(paste0(timing_path,blocks[1]), recursive = FALSE, full.names = FALSE)
  # for all blocks
  df_block = list()
  for (block in blocks) {
    # for all classes
    df_class = list()
    for (meth_class in meth_classes) {
      df_class[[meth_class]] <- do.call(rbind,lapply(list.files(paste0(timing_path, block, "/", meth_class),
                                                                pattern='REAL'),
                                                     function (time_file) {
                                                       fs <- strsplit(time_file, "_")[[1]][3]
                                                       data <- strsplit(time_file, "_")[[1]][1]
                                                       meth <- strsplit(strsplit(time_file, "_")[[1]][4], ".rds")[[1]][1]
                                                       time <- readRDS(paste0(timing_path, block, "/", meth_class, "/", time_file))
                                                       n_sample = ifelse(data=="BlREAL1",5,
                                                                         ifelse(data=="BlREAL2",12,
                                                                                ifelse(data=="BlREAL3",335,
                                                                                       ifelse(data=="SkREAL",4,
                                                                                              warning("You should specify the number of samples of your dataset.")))))
                                                       time = time/n_sample
                                                       if (is.null(time)) time <- NA
                                                       return(data.frame(values = time,
                                                                         score = "time",
                                                                         deconv = meth,
                                                                         dataset = data,
                                                                         feat_selec = fs,
                                                                         class=meth_class))
                                                     }))
    }
    df_block[[block]] = do.call(rbind,df_class)
    df_block[[block]]$block = block
  }
  saveRDS(do.call(rbind,df_block), paste0(score_path, "vivo_time.rds"))
}
