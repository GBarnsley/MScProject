#'
#'@export
MeanSquareJumpDistance <- function(matrixOfSamples){
  diffs <- diff(matrixOfSamples)
  for(i in 1:nrow(diffs)){
    diffs[i,] <- cumsum(diffs[i,])
  }
  diffs <- rowSums(diffs)^2
  return(mean(diffs))
}