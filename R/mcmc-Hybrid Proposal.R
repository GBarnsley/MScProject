#'
#' @export
hybridProp <- function(epiModel, hyperParameters, i){
  UseMethod("hybridProp")
}
#'
#' @export
hybridProp.SIR <- function(epiModel, hyperParameters, i){
  if(i%%3 == 0){
    randomWalkBeta(
      epiModel,
      hyperParameters
    )
  }
  else if(i%%3 == 1){
    randomWalkBeta(
      epiModel,
      hyperParameters
    )
  }
  else{
    stepProp(
      epiModel,
      hyperParameters,
      i
    )
  }
}
