# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
#
#       Options:
#       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii)  QN = Quantile normalization of input mixture (default = TRUE)
#       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
#               - note that cell subsets will be scaled by their absolute levels and will not be
#                 represented as fractions (to derive the default output, normalize absolute
#                 levels such that they sum to 1 for each mixture sample)
#               - the sum of all cell subsets in each mixture sample will be added to the output
#                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
#       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#               - sig.score = for each mixture sample, define S as the median expression
#                 level of all genes in the signature matrix divided by the median expression
#                 level of all genes in the mixture. Multiple cell subset fractions by S.
#               - no.sumto1 = remove sum to 1 constraint
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


#dependencies
library(e1071)
library(parallel)
library(preprocessCore)

#Core algorithm
CoreAlg <- function(X, y, absolute, abs_method){
    
    #try different values of nu
    svn_itor <- 3
    
    res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
    }
    
    if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
        out <- mclapply(1:svn_itor, res, mc.cores=6)
    
    nusvm <- rep(0,svn_itor)
    corrv <- rep(0,svn_itor)
    
    #do cibersort
    t <- 1
    while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- cor(k, y)
        t <- t + 1
    }
    
    #pick best model
    rmses <- nusvm
    mn <- which.min(rmses)
    model <- out[[mn]]
    
    #get and normalize coefficients
    q <- t(model$coefs) %*% model$SV
    q[which(q<0)]<-0
    if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
    if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
    
    mix_rmse <- rmses[mn]
    mix_r <- corrv[mn]
    
    newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
    
}

#do permutations
doPerm <- function(perm, X, Y, absolute, abs_method){
    itor <- 1
    Ylist <- as.list(data.matrix(Y))
    dist <- matrix()
    
    while(itor <= perm){
        #print(itor)
        
        #random mixture
        yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
        
        #standardize mixture
        yr <- (yr - mean(yr)) / sd(yr)
        
        #run CIBERSORT core algorithm
        result <- CoreAlg(X, yr, absolute, abs_method)
        
        mix_r <- result$mix_r
        
        #store correlation
        if(itor == 1) {dist <- mix_r}
        else {dist <- rbind(dist, mix_r)}
        
        itor <- itor + 1
    }
    newList <- list("dist" = dist)
}

#main function
CIBERSORT <- function(ref_profiles, dat, perm=0, QN=F, absolute=F, abs_method='sig.score'){
    if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score') stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
    
    #read in data
    X <- ref_profiles
    Y <- dat
    #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
    dups <- nrow(Y) - length(unique(rownames(Y)))
    if(dups > 0) {
        warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
        rownames(Y) <- make.names(Y[,1], unique=TRUE)
    }
    X <- data.matrix(X)
    Y <- data.matrix(Y)
    
    #order
    X <- X[order(rownames(X)),]
    Y <- Y[order(rownames(Y)),]
    
    #anti-log if max < 50 in mixture file
    if(max(Y) < 50) {Y <- 2^Y}
    
    #quantile normalization of mixture file
    if(QN == TRUE){
        tmpc <- colnames(Y)
        tmpr <- rownames(Y)
        Y <- normalize.quantiles(Y)
        colnames(Y) <- tmpc
        rownames(Y) <- tmpr
    }
    
    #store original mixtures
    Ymedian <- max(median(Y),1)
    
    #intersect genes
    Xgns <- row.names(X)
    Ygns <- row.names(Y)
    XYinter <- intersect(Xgns,Ygns)
    Y <- Y[XYinter,]
    X <- X[XYinter,]
    
    #standardize sig matrix
    X <- (X - mean(X))/sd(as.vector(X))
    
    #empirical null distribution of correlation coefficients
    if(perm > 0) {nulldist <- sort(doPerm(perm, X, Y, absolute, abs_method)$dist)}
    header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
    if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
    output <- matrix()
    itor <- 1
    samples <- ncol(Y)
    pval <- 9999
    
    #iterate through samples
    while(itor <= samples){
        #if (round(itor/30)==itor/30) {print(paste0("CIBERSORT progress ",round(100*itor/samples),"%"))}
        y <- Y[,itor]
        
        #standardize mixture
        y <- (y - mean(y))/sd(y)
        
        #run SVR core algorithm
        result <- CoreAlg(X, y, absolute, abs_method)
        
        #get results
        w <- result$w
        mix_r <- result$mix_r
        mix_rmse <- result$mix_rmse
        if(absolute && abs_method == 'sig.score') {
            w <- w * median(Y[,itor]) / Ymedian
        }
        
        #calculate p-value
        if(perm > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
        
        #print output
        out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
        if(absolute) out <- c(out, sum(w))
        if(itor == 1) {output <- out}
        else {output <- rbind(output, out)}
        itor <- itor + 1
    }
    
    #return matrix object containing all results
    obj <- rbind(header,output)
    obj <- obj[,-1]
    obj <- obj[-1,]
    obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
    rownames(obj) <- colnames(Y)
    if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
    else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
    return(obj)
}
