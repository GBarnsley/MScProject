#'
#'@export
acceptanceProb <- function(epiModel){
  UseMethod("acceptanceProb")
}
#'
#'@export
acceptanceProb.SIR <- function(epiModel){
  BetaProb <- sum(diff(epiModel@Samples[,'Beta'])!=0)/nrow(epiModel@Samples)
  GammaProb <- sum(diff(epiModel@Samples[,'Gamma'])!=0)/nrow(epiModel@Samples)
return(c(BetaProb, GammaProb))
}
#'
#'@export
acceptanceProb.iSIR <- function(epiModel){
  BetaGammaProbs <- acceptanceProb.SIR(epiModel)
  StepProb <- sum(rowSums(diff(epiModel@Samples[,-c(1,2)]) != 0)!=0)/(nrow(epiModel@Samples) - 1)
  return(c(BetaGammaProbs, StepProb))
}
acceptanceProb.rSIR <- function(epiModel){
  return(acceptanceProb.iSIR(epiModel))
}
