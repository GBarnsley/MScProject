#' Log-Likelihood
#'
#' Given a model with no unknown parameters this calculates the natural logarithm
#' of the likelihood. Only accepts SIR models.
#'
#' @param epiModel An epidemic model with know parameters
#' @return the log-likelihood
#' @export
logLikelihood <- function(epiModel){
  UseMethod("logLikelihood")
}
#' Log-likelihood for the SIR class
#' @export
logLikelihood.SIR <- function(epiModel){
  ll <- 0
  for(i in 1:length(epiModel@newI)){
    ll <- ll +
      dbinom(epiModel@newI[i],
             epiModel@S[i],
             probGen(
               epiModel@Beta*
                 epiModel@I[i]*
                 epiModel@t.step),
             log = TRUE) +
      dbinom(epiModel@newR[i],
             epiModel@I[i],
             probGen(
               epiModel@Gamma*
                 epiModel@t.step),
             log = TRUE)
    if(is.infinite(ll)){
      return(-Inf)
    }
  }
  return(ll)
}
