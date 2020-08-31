#' Function for calculating the mean squared jump distance from a matrix of
#' samples. Please see the written report for details of the method.
#' @param matrixOfSamples A matrix of samples
#' @return The mean squared jump distance
#' @export
MeanSquareJumpDistance <- function(matrixOfSamples){
  diffs <- diff(matrixOfSamples)
  for(i in 1:nrow(diffs)){
    diffs[i,] <- cumsum(diffs[i,])
  }
  diffs <- rowSums(diffs)^2
  return(mean(diffs))
}
