#'
#'@export
AcceptanceProbability <- function(samples){
  UseMethod("AcceptanceProbability")
}
#'
#'@export
AcceptanceProbability.numeric <- function(samples){
  return(mean(diff(samples) != 0))
}
#'
#'@export
AcceptanceProbability.matrix <- function(samples){
  return(mean(rowSums(diff(samples)^2) != 0))
}