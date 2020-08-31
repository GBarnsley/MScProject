#' Generic function for calculating acceptance probabilities from a set of samples
#' @param samples A matrix or vector of samples
#' @return The acceptance probability
#' @export
AcceptanceProbability <- function(samples){
  UseMethod("AcceptanceProbability")
}
#' Function for calculating acceptance probabilities from a vector of samples.
#' Checks if each element is different from is neighbour and then finds the mean
#' of that.
#' @param samples A vector of samples
#' @return The acceptance probability
#' @export
AcceptanceProbability.numeric <- function(samples){
  return(mean(diff(samples) != 0))
}
#' Function for calculating acceptance probabilities from a matrix of samples.
#' Checks if each row is different from is neighbour and then finds the mean
#' of that.
#' @param samples A matrix of samples
#' @return The acceptance probability
#' @export
AcceptanceProbability.matrix <- function(samples){
  return(mean(rowSums(diff(samples)^2) != 0))
}
