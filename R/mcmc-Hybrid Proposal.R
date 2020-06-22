#'
#' @export
hybridProp <- function(epiModel, hyperParameters){
  UseMethod("hybridProp")
}
#'
#' @export
hybridProp.SIR <- function(epiModel, hyperParameters){
  epiModel <- randomWalk.SIR(epiModel, hyperParameters)
  if(identical(epiModel,NA)){
    return(NA)
  }
  return(
    stepProp(epiModel, hyperParameters)
  )
}
#'
#' @export
logHybrid <- function(epiModel1, epiModel2){
  UseMethod("logHybrid")
}
#'
#' @export
logHybrid.SIR <- function(epiModel1, epiModel2){
  return(
    logStep(epiModel1, epiModel2)
  )
}
