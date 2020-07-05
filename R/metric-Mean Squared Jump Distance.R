#'
#'@export
meanSquaredJumping <- function(epiModel){
  distances <- diff(epiModel@Samples[,-c(1,2)])
  for(i in 1:nrow(distances)){
    distances[i,] <- cumsum(distances[i,])
  }
  distances <- rowSums(distances)
  return(mean(distances^2))
}
