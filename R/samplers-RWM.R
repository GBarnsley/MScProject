#'
#'@export
sampler_RW_block_TRACER_setup <- function(model, mvSaved, target, control) {
  #Tracer Columns
  evalCol <- as.integer(control$Column)
  ## control list extraction
  adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
  adaptScaleOnly      <- if(!is.null(control$adaptScaleOnly))      control$adaptScaleOnly      else FALSE
  adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 200
  adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
  scale               <- if(!is.null(control$scale))               control$scale               else 1
  propCov             <- if(!is.null(control$propCov))             control$propCov             else 'identity'
  tries               <- if(!is.null(control$tries))               control$tries               else 1
  ## node list generation
  targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
  calcNodes <- model$getDependencies(target)
  finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
  if(!is.integer(finalTargetIndex) |
     length(finalTargetIndex) != 1 |
     is.na(finalTargetIndex[1]))
    stop('Problem with target node in sampler_RW_block')
  calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
  calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
  ##calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
  isStochCalcNodesDepStage <- model$isStoch(calcNodesDepStage)   ## should be made faster
  calcNodesDepStageDeterm <- calcNodesDepStage[!isStochCalcNodesDepStage]
  calcNodesDepStageStoch <- calcNodesDepStage[isStochCalcNodesDepStage]
  ## numeric value generation
  scaleOriginal <- scale
  timesRan      <- 0
  timesAccepted <- 0
  timesAdapted  <- 0
  d <- length(targetAsScalar)
  scaleHistory  <- c(0, 0)                                                 ## scaleHistory
  acceptanceHistory  <- c(0, 0)                                            ## scaleHistory
  propCovHistory <- if(d<=10) array(0, c(2,d,d)) else array(0, c(2,2,2))   ## scaleHistory
  saveMCMChistory <- if(nimbleOptions('MCMCsaveHistory')) TRUE else FALSE
  if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
  propCovOriginal <- propCov
  chol_propCov <- chol(propCov)
  chol_propCov_scale <- scale * chol_propCov
  empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
  ## nested function and function list definitions
  ##my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
  targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
  ##my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)   ## old syntax: missing target argument
  my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
  ## checks
  if(!inherits(propCov, 'matrix'))        stop('propCov must be a matrix\n')
  if(!inherits(propCov[1,1], 'numeric'))  stop('propCov matrix must be numeric\n')
  if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
  if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
}
#'
#'@export
sampler_RW_block_TRACER_run <- function() {
  for(i in 1:tries) {
    currentValueVector <- values(model, targetNodesAsScalar)
    propValueVector <- generateProposalVector()
    ##lpMHR <- my_setAndCalculateDiff$run(propValueVector)
    values(model, targetNodesAsScalar) <<- propValueVector
    lpD <- calculateDiff(model, calcNodesProposalStage)
    if(lpD == -Inf) {
      nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      ## Drawing a random number is needed during first testing
      ## of this step in order to keep the random numbers identical
      ## to old behavior to see if tests that depend on particular
      ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
      ## we call decide() to ensure same behavior.
      ## jump <- decide(lpD)
      ## When new behavior is acceptable, we can remove the above line
      ## and uncomment the following:
      jump <- FALSE
    } else {
      ##jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## will use lpMHR - 0
      lpD <- lpD + calculateDiff(model, calcNodesDepStage)
      jump <- decide(lpD)
      if(jump) {
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage)#, logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageDeterm)#, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageStoch)#, logProbOnly = TRUE)
        model[["tracers"]][1:2,evalCol] <<- c(model[["tracers"]][1,evalCol] + 1,
                                              model[["tracers"]][2,evalCol] +
                                                sum((propValueVector - currentValueVector)^2))s
        
      } else {
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage)#, logProb = TRUE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesDepStageDeterm)#, logProb = FALSE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesDepStageStoch)#, logProbOnly = TRUE)
      }
    }
    if(adaptive)     adaptiveProcedure(jump)
  }
}
#'
#'@export
sampler_RW_block_TRACER_methods <- list(
  generateProposalVector = function() {
    propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
    returnType(double(1))
    return(propValueVector)
  },
  adaptiveProcedure = function(jump = logical()) {
    timesRan <<- timesRan + 1
    if(jump)     timesAccepted <<- timesAccepted + 1
    if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
    if(timesRan %% adaptInterval == 0) {
      acceptanceRate <- timesAccepted / timesRan
      timesAdapted <<- timesAdapted + 1
      if(saveMCMChistory) {
        setSize(scaleHistory, timesAdapted)                 ## scaleHistory
        scaleHistory[timesAdapted] <<- scale                ## scaleHistory
        setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
        acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        if(d <= 10) {
          propCovTemp <- propCovHistory                                           ## scaleHistory
          setSize(propCovHistory, timesAdapted, d, d)                             ## scaleHistory
          if(timesAdapted > 1)                                                    ## scaleHistory
            for(iTA in 1:(timesAdapted-1))                                      ## scaleHistory
              propCovHistory[iTA, 1:d, 1:d] <<- propCovTemp[iTA, 1:d, 1:d]    ## scaleHistory
          propCovHistory[timesAdapted, 1:d, 1:d] <<- propCov[1:d, 1:d]            ## scaleHistory
        }
      }
      adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
      scale <<- scale * adaptFactor
      ## calculate empirical covariance, and adapt proposal covariance
      if(!adaptScaleOnly) {
        gamma1 <- my_calcAdaptationFactor$getGamma1()
        for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
        empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
        propCov <<- propCov + gamma1 * (empirCov - propCov)
        chol_propCov <<- chol(propCov)
      }
      chol_propCov_scale <<- chol_propCov * scale
      timesRan <<- 0
      timesAccepted <<- 0
    }
  },
  getScaleHistory = function() {  ## scaleHistory
    if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
    returnType(double(1))
    return(scaleHistory)
  },
  getAcceptanceHistory = function() {  ## scaleHistory
    returnType(double(1))
    if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
    return(acceptanceHistory)
  },
  getPropCovHistory = function() { ## scaleHistory
    if(!saveMCMChistory | d > 10)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC and note that to reduce memory use we only save the proposal covariance history for parameter vectors of length 10 or less")
    returnType(double(3))
    return(propCovHistory)
  },
  ##getScaleHistoryExpanded = function() {                                                              ## scaleHistory
  ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)                         ## scaleHistory
  ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
  ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
  ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]                      ## scaleHistory
  ##    returnType(double(1)); return(scaleHistoryExpanded) },                                          ## scaleHistory
  ##getPropCovHistoryExpanded = function() {                                                            ## scaleHistory
  ##    propCovHistoryExpanded <- array(dim=c(timesAdapted*adaptInterval,d,d), init=FALSE)              ## scaleHistory
  ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
  ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
  ##            propCovHistoryExpanded[(iTA-1)*adaptInterval+j,1:d,1:d] <- propCovHistory[iTA,1:d,1:d]  ## scaleHistory
  ##    returnType(double(3)); return(propCovHistoryExpanded) },                                        ## scaleHistory
  reset = function() {
    scale   <<- scaleOriginal
    propCov <<- propCovOriginal
    chol_propCov <<- chol(propCov)
    chol_propCov_scale <<- chol_propCov * scale
    timesRan      <<- 0
    timesAccepted <<- 0
    timesAdapted  <<- 0
    if(saveMCMChistory) {
      scaleHistory  <<- c(0, 0)    ## scaleHistory
      acceptanceHistory  <<- c(0, 0)
      if(d <= 10)
        propCovHistory <<- nimArray(0, dim = c(2,d,d))
    }
    my_calcAdaptationFactor$reset()
  }
)