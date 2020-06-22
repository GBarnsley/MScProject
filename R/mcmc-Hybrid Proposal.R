#'
#' @export
hybridProp <- function(epiModel, hyperParameters, i){
  UseMethod("hybridProp")
}
#'
#' @export
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
}
