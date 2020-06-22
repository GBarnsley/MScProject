#' function that sets up initial values for MH algorithm.
#' @export
initialValues <- function(epiModel, hyperParameters){
  UseMethod("initialValues", epiModel)
}
#' Intial values for SIR model. Just sets Beta and Gamma to
#' the start value specified in the hyperparameters.
#' from it.
#' @export
initialValues.SIR <- function(epiModel, hyperParameters){
  epiModel@Beta <- hyperParameters$Beta
  epiModel@Gamma <- hyperParameters$Gamma
  return(
    epiModel
  )
}
#' Intial values for SIR model. Just sets Beta and Gamma to
#' the start value specified in the hyperparameters and then
#' generates newI.
#' from it.
#' @export
initialValues.iSIR <- function(epiModel, hyperParameters){
  epiModel@Beta <- hyperParameters$Beta
  epiModel@Gamma <- hyperParameters$Gamma
  return(
    simulate(epiModel)
  )
}
#' Intial values for SIR model. Just sets Beta and Gamma to
#' the start value specified in the hyperparameters and then
#' generates newR.
#' from it.
#' @export
initialValues.rSIR <- function(epiModel, hyperParameters){
  epiModel@Beta <- hyperParameters$Beta
  epiModel@Gamma <- hyperParameters$Gamma
  return(
    simulate(epiModel)
  )
}
