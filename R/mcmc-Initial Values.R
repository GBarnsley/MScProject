#' Generic function that calls the model specific initialization methods
#' @param epiModel An object of the class one of the specific epidemic model
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return An object of the given epidemic class but initialized via the specified method
#' @export
initialValues <- function(epiModel, hyperParameters){
  sampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = stepSampler_setup,
    run = stepSampler_run,
    methods = list(
      reset = function() {}
    )
  )
  UseMethod("initialValues", epiModel)
}
