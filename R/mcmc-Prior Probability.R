#' Calculates the log of the prior probability for the model
#' @param epiModel some epidemic model with no missing data
#' @param priors the hyperparameters of the priors in list form
#' @return the value of the log prior probability
#' @export
logPrior <- function(epiModel, priors){
  UseMethod("logPrior")
}
#' Log prior for when newI is a parameter in the SIR model.
#' The priors from Beta and Gamma are considered to be Gamma(alpha, beta).
#' The prior for newI(t) is the discrete Unif(max(newR(t+1) - I(t), 0), S(t)).
#' @export
logPrior.iSIR <- function(epiModel, priors){
  output <-
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
  #for(t in 1:length(epiModel@S[epiModel@S!=0])){
  #  output <-
  #    output -
  #    log(
  #      epiModel@S[t] -
  #        max(
  #          epiModel@newR[t+1] -
  #            epiModel@I[t],
  #          0,
  #          na.rm = TRUE
  #        )
  #    )
  #}
  return(
    output
  )
}
#' Log prior for when newR is a parameter in the SIR model.
#' The priors from Beta and Gamma are considered to be Gamma(alpha, beta).
#' The prior for newR(t) is the discrete Unif(0, I(t)).
#' Note, I am sure the prior for newR(t) can be more infromative, e.g.
#' needing atleast one infectious person until the end etc. I've not
#' gotten round to thinking about it yet, if I ever do.
#' @export
logPrior.rSIR <- function(epiModel, priors){
  output <-
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
    )# -
    #sum(log(epiModel@I[epiModel@I!=0]))
  return(
    output
  )
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
