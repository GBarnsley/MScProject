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
  for(i in 1:burnin){
    current <- proposal(current, hyperParameters, i)
  }
  count <- samples
  for(i in 1:samples){
    for(j in 1:keep){
      current <- proposal(current, hyperParameters, j + (i-1)*keep)
    }
    result[[i]] <- current
    cat(count,"\r")
    count <- count-1
  }
  cat("\r", "Probability of Acceptance: ", result[[samples]]@acceptance/result[[samples]]@attempts*100, "%\n")
  class(result) <- c("samples", class(epiModel))
  return(
    result
  )
}
