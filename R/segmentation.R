# file:   segmentation.R
# author: Gavin Ha, Ph.D.
#         Fred Hutchinson Cancer Research Center
# contact: <gha@fredhutch.org>
# website: https://GavinHaLab.org
#
# ichorCNA website: https://github.com/GavinHaLab/ichorCNA
# date:   January 6, 2020
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.
#' @export
HMMsegment <- function(x, validInd = NULL, dataType = "copy", param = NULL, 
    chrTrain = c(1:22), maxiter = 50, estimateNormal = TRUE, estimatePloidy = TRUE, 
    estimatePrecision = TRUE, estimateVar = TRUE, estimateSubclone = TRUE, estimateTransition = TRUE,
    estimateInitDist = TRUE, logTransform = FALSE, likChangeConvergence = 1e-3, verbose = TRUE) {
  	chr <- as.factor(seqnames(x[[1]]))
	# setup columns for multiple samples #
	dataMat <- as.matrix(as.data.frame(lapply(x, function(y) { mcols(y)[, dataType] })))
	
	# normalize by median and log data #
	if (logTransform){
    dataMat <- apply(dataMat, 2, function(x){ log(x / median(x, na.rm = TRUE)) })
	}else{
	  dataMat <- log(2^dataMat)
	}
	## update variable x with loge instead of log2
  for (i in 1:length(x)){
    mcols(x[[i]])[, dataType] <- dataMat[, i]
  }
  if (!is.null(chrTrain)) {
		chrInd <- chr %in% chrTrain
  }else{
  	chrInd <- !logical(length(chr))
  }
  if (!is.null(validInd)){
    chrInd <- chrInd & validInd
  }  

	if (is.null(param)){
		param <- getDefaultParameters(dataMat[chrInd])
	}
	#if (param$n_0 == 0){
	#	param$n_0 <- .Machine$double.eps
	#}
	####### RUN EM ##########
  convergedParams <- runEM(dataMat, chr, chrInd, param, maxiter, 
      verbose, estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, 
      estimateVar = estimateVar, estimateSubclone = estimateSubclone, 
      estimatePrecision = estimatePrecision, estimateTransition = estimateTransition, 
      estimateInitDist = estimateInitDist, likChangeConvergence = likChangeConvergence)
	
  # Calculate likelihood using converged params
 # S <- param$numberSamples
 # K <- length(param$ct)
 # KS <- K ^ S
 # py <- matrix(0, KS, nrow(dataMat))
 # iter <- convergedParams$iter
  # lambdasKS <- as.matrix(expand.grid(as.data.frame(convergedParams$lambda[, , iter])))
  # for (ks in 1:KS) {
  #   probs <- tdistPDF(dataMat, convergedParams$mus[ks, , iter], lambdasKS[ks, ], param$nu)
  #   py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
  # }
  # 
  viterbiResults <- runViterbi(convergedParams, chr)
  
  # setup columns for multiple samples #
  segs <- segmentData(x, validInd, viterbiResults$states, convergedParams)
  #output$segs <- processSegments(output$segs, chr, start(x), end(x), x$DataToUse)
  names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0("HLAMP", 2:1000))
  #if (c(0) %in% param$ct){ #if state 0 HOMD is IN params#
  	#names <- c("HOMD", names)
  	# shift states to start at 2 (HETD)
    #tmp <- lapply(segs, function(x){ x$state <- x$state + 1; x})
    #viterbiResults$states <- as.numeric(viterbiResults$states) + 1
	#}
	### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
  cnaList <- list()
  S <- length(x)
  for (s in 1:S){
    id <- names(x)[s]
    copyNumber <- param$jointCNstates[viterbiResults$state, s]
    subclone.status <- param$jointSCstatus[viterbiResults$state, s]
  	cnaList[[id]] <- data.frame(cbind(sample = as.character(id), 
                  chr = as.character(seqnames(x[[s]])),	
                  start = start(x[[s]]), end = end(x[[s]]), 
                  copy.number = copyNumber,
                  event = names[copyNumber + 1], 
                  logR = round(log2(exp(dataMat[,s])), digits = 4),
                  subclone.status = as.numeric(subclone.status)
  	))
  
    cnaList[[id]] <- transform(cnaList[[id]], 
                              start = as.integer(as.character(start)),
                              end = as.integer(as.character(end)), 
                              copy.number = as.numeric(copy.number),
                              logR = as.numeric(as.character(logR)),
                              subclone.status = as.numeric(subclone.status))
  
  	## order by chromosome ##
  	chrOrder <- unique(chr) #c(1:22,"X","Y")
  	cnaList[[id]] <- cnaList[[id]][order(match(cnaList[[id]][, "chr"],chrOrder)),]
  	## remove MT chr ##
    cnaList[[id]] <- cnaList[[id]][cnaList[[id]][,"chr"] %in% chrOrder, ]
    
    ## segment mean loge -> log2
    #segs[[s]]$median.logR <- log2(exp(segs[[s]]$median.logR))
    segs[[s]]$median <- log2(exp(segs[[s]]$median))
    ## add subclone status
    segs[[s]]$subclone.status <-  param$jointSCstatus[segs[[s]]$state, s]
  }	
  convergedParams$segs <- segs
  return(list(cna = cnaList, results = convergedParams, viterbiResults = viterbiResults))
}
#' @export
getTransitionMatrix <- function(K, e, strength){
  A <- matrix(0, K, K)
  for (j in 1:K) {
    A[j, ] <- (1 - e[1]) / (K - 1)
    A[j, j] <- e[1]
  }
  A <- normalize(A)
  A_prior <- A
  dirPrior <- A * strength[1]
  return(list(A=A, dirPrior=dirPrior))
}
#' @export
getDefaultParameters <- function(x, maxCN = 5, ct.sc = c(1,3), n_0 = 0.5, ploidy_0 = 2, 
                                 normal2IgnoreSC = 0.9, e = 0.9999999, e.subclone = 0.1,
                                 e.sameState = 50, strength = 10000000, 
                                 includeHOMD = FALSE, likModel = "t"){
  if (includeHOMD){
    ct <- 0:maxCN
  }else{
    ct <- 1:maxCN
  }
  #ct <- c(ct, maxCN * 5) ## add additional CN state to catch extreme HLAMP
	param <- list(
		strength = strength, e = e,
		ct = c(ct, ct.sc),
		ct.sc.status = c(rep(FALSE, length(ct)), rep(TRUE, length(ct.sc))),
		phi_0 = ploidy_0, alphaPhi = 4, betaPhi = 0.75,
		n_0 = n_0, alphaN = 1, betaN = 1,
		sp_0 = 0.5, alphaSp = 2, betaSp = 2,
		lambda = as.matrix(rep(100, length(ct)+length(ct.sc)), ncol=1),
		nu = 2.1,
		kappa = rep(75, length(ct)), 
		alphaLambda = 2,
		betaVar = 2,
		likModel = likModel
	)
	K <- length(param$ct)
	S <- ncol(x)
	KS <- K ^ S
	
	if (grepl("gauss", likModel, ignore.case = TRUE)){
	  likModel <- "Gaussian"
	  param$likModel <- likModel
	}
  ## initialize hyperparameters for precision using observed data ##
	# multiple or single samples (columns)
  param$numberSamples <- S
  #betaLambdaVal <- ((apply(x, 2, function(x){ sd(diff(x), na.rm=TRUE) }) / sqrt(length(param$ct))) ^ 2)
  
  logR.var <- ( apply(x, 2, sd, na.rm = TRUE) / sqrt(K) ) ^ 2
  #logR.var <- ( apply(x, 2, var, na.rm = TRUE) / sqrt(K) ) 
	
	#param$betaLambda <- betaLambdaVal #t(replicate(K, betaLambdaVal))
	param$betaLambda <- matrix(logR.var, ncol = param$numberSamples, nrow = K, byrow = TRUE)
  param$alphaLambda <- rep(param$alphaLambda, K)
  
	#logR.var <- 1 / ((apply(x, 2, sd, na.rm = TRUE) / sqrt(K)) ^ 2)
	logR.lambda <- 1 / logR.var
  if (ncol(x) > 1){ # multiple samples (columns)
		param$lambda <- matrix(logR.lambda, nrow=K, ncol=S, byrow=T, dimnames=list(c(),colnames(x)))
	}else{ # only 1 sample and using student's-t likelihood model
		#logR.var <- 1 / ((sd(x, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
    param$lambda <- matrix(logR.lambda, length(param$ct))
    param$lambda[param$ct %in% c(2)] <- logR.lambda * 10
    param$lambda[param$ct %in% c(1,3)] <- logR.lambda 
    param$lambda[param$ct >= 4] <- logR.lambda
    param$lambda[param$ct == max(param$ct)] <- logR.lambda / 15
    param$lambda[param$ct.sc.status] <- logR.lambda * 2
    param$betaLambda <- param$alphaLambda / param$lambda
    #param$betaLambda[param$ct %in% c(2)] <- param$betaLambda[param$ct %in% c(2)] 
    #param$betaLambda[param$ct.sc.status] <- param$betaLambda[param$ct.sc.status] 
  }

  ## initialize hyperparameters and variance parameters for Gaussians
  if (likModel == "Gaussian"){
    param$psi <- diag(param$numberSamples) # parameter for inverse-wishart prior
    param$nu
    if (numSamples > 1){
      covar <- getCovarianceMatrix(x)
      param$cor <- covar$cor
      param$covar <- covar$covar
      param$sw <- rep(1, param$numberSamples)
    }else{
      param$cor <- var(x)
      param$covar <- 1
      param$sw <- 1
    }
    #param$var <- t(replicate(K, diag(covar$covar)))
    #param$var <- matrix(apply(x, 2, var, na.rm=T), ncol = S, nrow = K, byrow = TRUE)
    param$betaVar <- rep(param$betaVar, S)
    #alphaVar <- 1 / (apply(x, 2, var, na.rm = TRUE) / sqrt(K) ^ 2)
    #alphaVar <- 1 / logR.var  # same as betaLambda
    param$alphaVar <- matrix(1 / logR.var, ncol = S, nrow = K, byrow = TRUE)
    param$var <- matrix(logR.var, ncol = S, nrow = K, byrow = TRUE)
    param$var[param$ct %in% c(2)] <- logR.var / 10
    param$var[param$ct %in% c(1,3)] <- logR.var 
    param$var[param$ct >= 4] <- logR.var
    param$var[param$ct == max(param$ct)] <- logR.var * 15
    param$var[param$ct.sc.status] <- logR.var * 2
    param$alphaVar <- param$betaVar / param$var
    #param$alphaVar[param$ct.sc.status] <- param$alphaVar[param$ct.sc.status] * 2

    ### custom settings
    # highest CN state has higher variance to capture outliers
    #param$alphaVar[param$jointCNstates == max(param$ct)] <- param$alphaVar[param$jointCNstates == max(param$ct)] * 10
    # subclonal states ahve higher variance to capture more subclonal range
    #param$alphaVar[param$ct.sc.status, 1] <- param$alphaVar[param$ct.sc.status, 1] * 2
  }

  # define joint copy number states #
	param$jointStates <- expand.grid(rep(list(1:length(param$ct)), S))
  param$jointCNstates <- expand.grid(rep(list(param$ct), S))
  param$jointSCstatus <- expand.grid(rep(list(param$ct.sc.status), S))
  #KS <- nrow(param$jointCNstates)
  #param$covar <- replicate(KS, covar$covar)
	
  colnames(param$jointCNstates) <- paste0("Sample.", 1:param$numberSamples)
  colnames(param$jointSCstatus) <- paste0("Sample.", 1:param$numberSamples)
  
	# Initialize transition matrix to the prior
	txn <- getTransitionMatrix(K ^ S, e, strength)
  ## set higher transition probs for same CN states across samples ##
  # joint states where at least "e.sameState" fraction of samples with the same CN state
	#apply(param$jointCNstates, 1, function(x){ sum(duplicated(as.numeric(x))) > 0 })
  cnStateDiff <- apply(param$jointCNstates, 1, function(x){ (abs(max(x) - min(x)))})
  if (e.sameState > 0 & S > 1){
		txn$A[, cnStateDiff == 0] <- txn$A[, cnStateDiff == 0] * e.sameState * K 
		#txn$A[, cnStateDiff >= 3] <- txn$A[, cnStateDiff >=3]  / e.sameState / K
	}
  param$kappa <- rep(75, K ^ S)
  param$kappa[cnStateDiff == 0] <- param$kappa[cnStateDiff == 0] + 125
  param$kappa[cnStateDiff >=3] <- param$kappa[cnStateDiff >=3] - 50
  param$kappa[which(rowSums(param$jointCNstates==2) == S)] <- 800
  # penalize homozygous deletion state
  if (includeHOMD){
    K <- length(param$ct)
    txn$A[1, 2:K] <- txn$A[1, 2:K] * 1e-5; 
    txn$A[2:K, 1] <- txn$A[2:K, 1] * 1e-5;
    txn$A[1, 1] <- txn$A[1, 1] * 1e-5
  }
  # penalize transitions into subclonal states
  #ind.bothSC <- rowSums(param$jointSCstatus) == S
  if (e.subclone > 1){
    ind.SC <- rowSums(param$jointSCstatus) > 0
    txn$A[, ind.SC] <- txn$A[, ind.SC] / (e.subclone * K)
    # penalize transitions into outlier maxCN state
    ind.maxCN <- rowSums(param$jointCNstates == max(param$ct)) > 0
    #txn$A[, ind.maxCN] <- txn$A[, ind.maxCN] / (2 * K)
  }
  txn$A <- normalize(txn$A)
	param$A <- txn$A
  param$A_prior <- param$A
  param$dirPrior <- param$A * strength[1] 

  # ignore subclones for sample with n >= normal2IgnoreSC
  indKS <- 1:nrow(param$jointSCstatus) # use all joint states
  indK <- 1:length(param$ct.sc.status)
  if (sum(n_0 >= normal2IgnoreSC) > 0){
    for (i in which(n_0 >= normal2IgnoreSC)){
      # keep states excluding subclones for samples with n >= normal2IgnoreSC
      indKS <- intersect(indKS, which(!param$jointSCstatus[, i]))
      indK <- intersect(indK, which(param$ct.sc.status))
      param$lambda[indK, i] <- NaN
      param$betaLambda[indK, i] <- NaN
      if (likModel == "Gaussian"){
        param$var[indK, i] <- NaN
      }    
    }

    param$kappa <- param$kappa[indKS]
    param$alphaVar <- param$alphaVar[indK, , drop = FALSE]
    param$jointStates <- param$jointStates[indKS, , drop = FALSE]
    param$jointCNstates <- param$jointCNstates[indKS, , drop = FALSE]
    param$jointSCstatus <- param$jointSCstatus[indKS, , drop = FALSE]
    param$A <- normalize(param$A[indKS, indKS])
    param$A_prior <- param$A
    param$dirPrior <- param$A * strength[1]
  } 
  
  return(param)
}


#' @export
segmentData <- function(dataGR, validInd, states, convergedParams){
  if (sum(convergedParams$param$ct == 0) ==0){
  	includeHOMD <- FALSE
  }else{
  	includeHOMD <- TRUE
  }
  if (!includeHOMD){
    names <- c("HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }else{
    names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }
  states <- states[validInd]
  S <- length(dataGR)
  jointStates <- convergedParams$param$jointCNstates
  jointSCstatus <- convergedParams$param$jointSCstatus
  colNames <- c("seqnames", "start", "end", "copy")
  segList <- list()
  for (i in 1:S){
  	id <- names(dataGR)[i]
  	dataIn <- dataGR[[i]][validInd, ]
    rleResults <- t(sapply(runValue(seqnames(dataIn)), function(x){
      ind <- as.character(seqnames(dataIn)) == x
      r <- rle(states[ind])
    }))
    rleLengths <- unlist(rleResults[, "lengths"])
    rleValues <- unlist(rleResults[, "values"])
    sampleDF <- as.data.frame(dataIn)
    numSegs <- length(rleLengths)
    segs <- as.data.frame(matrix(NA, ncol = 7, nrow = numSegs, 
                   dimnames = list(c(), c("chr", "start", "end", "state", "event", "median", "copy.number"))))
    prevInd <- 0
    for (j in 1:numSegs){
      start <- prevInd + 1
      end <- prevInd + rleLengths[j]
      segDF <- sampleDF[start:end, colNames]
      prevInd <- end
      numR <- nrow(segDF)
      segs[j, "chr"] <- as.character(segDF[1, "seqnames"])
      segs[j, "start"] <- segDF[1, "start"]
      segs[j, "state"] <- rleValues[j]
      segs[j, "copy.number"] <- jointStates[rleValues[j], i]
      if (segDF[1, "seqnames"] == segDF[numR, "seqnames"]){
        segs[j, "end"] <- segDF[numR, "end"]
        segs[j, "median"] <- round(median(segDF$copy, na.rm = TRUE), digits = 6)
        if (includeHOMD){
        	segs[j, "event"] <- names[segs[j, "copy.number"] + 1]
        }else{
        	segs[j, "event"] <- names[segs[j, "copy.number"]]
        }
      }else{ # segDF contains 2 different chromosomes
        print(j)
      }                                      
    }
    segList[[id]] <- segs
  }
  return(segList)
}
    
#' @export
runViterbi <- function(convergedParams, chr){
  message("runViterbi: Segmenting and classifying")
  chrs <- levels(chr)
  chrsI <- vector('list', length(chrs))
  # initialise the chromosome index and the init state distributions
  for(i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }
  segs <- vector('list', length(chrs))
  py <- convergedParams$py
  N <- ncol(py)
  Z <- rep(0, N)
  convergeIter <- convergedParams$iter
  piG <- convergedParams$pi[, convergeIter]
  A <- convergedParams$A


  for(c in 1:length(chrsI)) {
    I <- chrsI[[c]]
    output <- .Call("viterbi", log(piG), log(A), log(py[, I]), PACKAGE = "HMMcopy")
    Z[I] <- output$path
    segs[[c]] <- output$seg
  }
  return(list(segs=segs, states=Z))
}
#' @export
# Normalize a given array to sum to 1
normalize <- function(A) {
	vectorNormalize <- function(x){ x / (sum(x) + (sum(x) == 0)) }
	if (length(dim(A)) < 2){
  	M <- vectorNormalize(A)
  }else{
  	M <- t(apply(A, 1, vectorNormalize))
  }
  return(M);
}


# processSegments <- function(seg, chr, start, end, copy) {
#   segment <- data.frame()
#   chromosomes <- levels(chr)
#   for (i in 1:length(chromosomes)) {
#     seg_length = dim(seg[[i]])[1]
#     chr_name <- rep(chromosomes[i], seg_length)
#     chr_index <- which(chr == chromosomes[i])
#     chr_start <- start[chr_index][seg[[i]][, 1]]
#     chr_stop <- end[chr_index][seg[[i]][, 2]]
#     chr_state <- seg[[i]][, 3]
#     chr_median <- rep(0, seg_length)
#     for(j in 1:seg_length) {
#       chr_median[j] <-
#         median(na.rm = TRUE, log2(exp(copy[chr_index][seg[[i]][j, 1]:seg[[i]][j, 2]])))
#     }
#     segment <- rbind(segment, cbind(chr = chr_name,
#       start = as.numeric(chr_start), end = chr_stop, state = chr_state,
#       median = chr_median))
#   }
#   segment <- transform(segment, start = as.numeric(as.character(start)),
#     end = as.numeric(as.character(end)), as.numeric(as.character(state)),
#     median = as.numeric(as.character(median)))
#   return(segment)
# }
