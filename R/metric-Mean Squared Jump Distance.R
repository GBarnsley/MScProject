#'
#'@export
meanSquaredJumping <- function(epiModel){
  distances <- diff(epiModel@Samples[,-c(1,2)])
  for(i in 1:nrow(distances)){
    distances[i,] <- cumsum(distances[i,]) != 0
  }
  distances <- rowSums(distances)
  return(mean(distances^2))
}
