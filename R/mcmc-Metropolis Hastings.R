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
                               thin = 10){
  epiModel <- initialValues(epiModel, hyperParameters)
  epiModel@MCMC <- buildMCMCInternal(epiModel, hyperParameters)
  epiModel@MCMC$run(niter = samples*thin + burnin, thin = thin, nburnin = burnin)
  epiModel@Samples <- as.matrix(epiModel@MCMC$mvSamples)
  epiModel@Metrics <- epiModel@Model$tracers/(samples*thin + burnin)
  row.names(epiModel@Metrics) <- c("Acceptance Probability", "Mean Squared Jump Distance", "Mean Squared Error")
  if(identical(NULL, hyperParameters$`True Values`)){
    epiModel@Metrics[3,1:ncol(epiModel@Metrics)] <- rep(NA, ncol(epiModel@Metrics))
  }
  else{
    epiModel@Metrics[3,1:ncol(epiModel@Metrics)] <- sqrt(epiModel@Metrics[3,1:ncol(epiModel@Metrics)])
  }
  epiModel@Metrics[2,1:ncol(epiModel@Metrics)] <- sqrt(epiModel@Metrics[2,1:ncol(epiModel@Metrics)])
  return(epiModel)
}

