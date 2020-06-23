#'
#' @export
<<<<<<< HEAD
hybridProp <- function(epiModel, hyperParameters, i){
=======
hybridProp <- function(epiModel, hyperParameters){
>>>>>>> 1c714e3e41b618fd584d43c4f5ae32a0417f471a
  UseMethod("hybridProp")
}
#'
#' @export
<<<<<<< HEAD
hybridProp.SIR <- function(epiModel, hyperParameters, i){
  return(stepProp(
    randomWalkGamma(
      randomWalkBeta(
        epiModel,
        hyperParameters),
      hyperParameters),
    hyperParameters,
    i)
    )
=======
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
>>>>>>> 1c714e3e41b618fd584d43c4f5ae32a0417f471a
}
