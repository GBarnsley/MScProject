#' Generates the misisng values/estimates parameters using the
#' Metropolis Hastings algorithm. Returns a list of samples from
#' the posterior distribution. Currently only runs the hybrid algorithm
#' @param epiModel Some object with parameters/missing data to impute
#' @param proposal the proposal function
#' @param hyperParameters the hyperparameters for the priors, in list form + variance for randomwalk
#' @param samples the target number of samples to be left with
#' @param burnin the amount of samples to be dropped as burnin
#' @param keep the algorithm will keep 1 in keep of the generated samples
#' @export
metropolisHastings <- function(epiModel,
                               hyperParameters,
                               samples = 1000,
                               burnin = 500,
                               thin = 10,
                               Evaluate = FALSE,
                               trueValues = NA){
  epiModel <- initialValues(epiModel, hyperParameters)
  epiModel@MCMC <- buildMCMCInternal(epiModel, hyperParameters)
  if(Evaluate){
    epiModel@MCMC$run(niter = samples*thin + burnin)
    epiModel@Samples <- as.matrix(epiModel@MCMC$mvSamples)
    epiModel@Metrics <- list(
      `Acceptance Probabilities` = acceptanceProb(epiModel),
      `Mean Squared Jumping Distance` = meanSquaredJumping(epiModel)
    )
    if(!is.na(trueValues)){
      epiModel@Metrics$`Mean Mean Squared Error` <- meanSquaredError(epiModel,
                                                                     trueValues)
    }
    epiModel@Samples <- epiModel@Samples[(burnin+1):(samples*thin + burnin),]
    epiModel@Samples <- epiModel@Samples[seq(1, samples*thin, thin),]
    return(epiModel)
  }
  else{
    epiModel@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
    epiModel@Samples <- as.matrix(epiModel@MCMC$mvSamples)
    return(epiModel)
    }
}

