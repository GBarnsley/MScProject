#' Calculates the log of the prior probability for the model
#' @param epiModel some epidemic model with no missing data
#' @param priors the hyperparameters of the priors in list form
#' @return the value of the log prior probability
#' @export
logPrior <- function(epiModel, priors){
  UseMethod("logPrior")
}
#' Log prior for when Beta and Gamma are the paramters considered.
#' @export
logPrior.SIR <- function(epiModel, priors){
  return(
    dgamma(
      epiModel@Beta,
      priors$Beta$Alpha,
      priors$Beta$Beta,
      log = TRUE
    ) +
      dgamma(
        epiModel@Gamma,
        priors$Gamma$Alpha,
        priors$Gamma$Beta,
        log = TRUE
      )
  )
}
