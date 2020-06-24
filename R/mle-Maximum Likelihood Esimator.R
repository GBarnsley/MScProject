#' Maximum Likelihood Estimate
#'
#' Uses optim to make an MLE for the parameters of given epidemic model.
#' Uses L-BFGS-B so bounds can be provided as a list of vectors.
#' Currently accepts only SIR class models
#'
#' @param epiModel An epidemic model
#' @param bounds list of vectors of lower and upper bounds for each parameter to be estimated
#' @return the MLE
#' @export
MLE <- function(epiModel, bounds){
  UseMethod(
    "MLE"
    )
}
#' MLE for SIR class
#' @export
MLE.SIR <- function(epiModel, bounds = list(Beta = c(10^-10,0.1),
                                            Gamma = c(10^-10,0.1))){
  if(is.null(epiModel@newI)|is.null(epiModel@newR)){
    print("Error: Insufficient Data")
    return(NA)
  }
  func <- function(param, epiModel){
    epiModel@Beta = param[1]
    epiModel@Gamma  = param[2]
    res <- logLikelihood(epiModel)
    return(res)
  }
  return(
    optim(
      par = c(sum(bounds$Beta)/2,sum(bounds$Gamma)/2),
      fn = func,
      epiModel = epiModel,
      method = "L-BFGS-B",
      lower = c(bounds$Beta[1],bounds$Gamma[1]),
      upper = c(bounds$Beta[2],bounds$Gamma[2]),
      control = list(
        fnscale = -1,
        trace = 0
      )
    )$par
  )
}
