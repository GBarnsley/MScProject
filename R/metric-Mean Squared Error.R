#'
#'@export
meanSquaredError <- function(epiModel, trueValues){
  UseMethod("meanSquaredError")
}
#'
#'@export
meanSquaredError.SIR <- function(epiModel, trueValues){
  Error <- c(mean(((epiModel@Samples[,'Beta'] - trueValues$Beta)^2)),
              mean((epiModel@Samples[,'Gamma'] - trueValues$Gamma)^2))
  return(Error)
}
#'
#'@export
meanSquaredError.iSIR <- function(epiModel, trueValues){
  Error <- meanSquaredError.SIR(epiModel, trueValues)
  newIError <- mean(rowMeans((epiModel@Samples[,-c(1,2)] - matrix(trueValues$newI,
                                      nrow = nrow(epiModel@Samples),
                                      ncol = length(trueValues$newI),
                                      byrow = TRUE))^2))
  return(c(Error, newIError))
}
#'
#'@export
meanSquaredError.rSIR <- function(epiModel, trueValues){
  return(meanSquaredError.iSIR(epiModel, trueValues))
}
