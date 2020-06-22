#' Generates the misisng values/estimates parameters using the
#' Metropolis Hastings algorithm. Returns a list of samples from
#' the posterior distribution. Currently accepts iSIR or rSIR objects
#' and proposal takes randomWalk or independenceSampler
#' @param epiModel Some object with parameters/missing data to impute
#' @param proposal the proposal function
#' @param hyperParameters the hyperparameters for the priors, in list form + variance for randomwalk
#' @param samples the target number of samples to be left with
#' @param burnin the amount of samples to be dropped as burnin
#' @param keep the algorithm will keep 1 in keep of the generated samples
#' @export
metropolisHastings <- function(epiModel,
                               proposal,
                               hyperParameters,
                               samples = 1000,
                               burnin = 1000,
                               keep = 100){
  result <- list()
  current <- NA
  while(identical(current, NA)){
    current <- initialValues(epiModel, hyperParameters$`Initial Values`)
  }
  count <- samples
  accepted <- 0
  for(i in 1:(samples*keep + burnin)){
    candidate <- proposal(current, hyperParameters)
    if(identical(candidate,NA)){
      acceptance <- -Inf
    }
    else{
      acceptance <-
        logPrior(candidate, hyperParameters$Priors) +
        logLikelihood(candidate) +
        logProposal(candidate, current, proposal, hyperParameters) -
        logPrior(current, hyperParameters$Priors) -
        logLikelihood(current) -
        logProposal(current, candidate, proposal, hyperParameters)
      #print(c(candidate@Beta, candidate@Gamma, logLikelihood(candidate)))
      #print(c(current@Beta, current@Gamma, logLikelihood(current)))
      #print(exp(acceptance))
      #print(c(candidate@Beta, candidate@Gamma))
      #print(c(logPrior(candidate, hyperParameters$Priors),
      #        logLikelihood(candidate),
      #        logProposal(candidate, current, proposal, hyperParameters),
      #        logPrior(current, hyperParameters$Priors),
      #        logLikelihood(current),
      #        logProposal(current, candidate, proposal, hyperParameters)))
    }
    if(log(runif(1)) < acceptance){
      current <- candidate
      if(i > burnin){
      accepted <- accepted + 1
      }
    }
    if(i > burnin & ((i-burnin)%%keep)==0){
      result[[(i-burnin)/keep]] <- current
      cat(count,"\r")
      count <- count-1
    }
  }
  cat("\r", "Probability of Acceptance: ", accepted/(samples*keep), "\n")
  class(result) <- c("samples", class(epiModel))
  return(
    result
  )
}
