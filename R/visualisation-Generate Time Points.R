#'
#'@export
generateTimepoints <- function(samples){
  UseMethod("generateTimepoints", samples[[1]])
}
#'
#'@export
generateTimepoints.iSIR <- function(samples){
  usefulNewILength <- min(which(samples[[round(length(samples)/2)]]@S==0)) - 1
  timePoints <- c(sample(1:round(usefulNewILength/3), 1),
                  sample(round(usefulNewILength/3):round(usefulNewILength*2/3), 1),
                  sample(round(usefulNewILength*2/3):usefulNewILength, 1))
  return(timePoints)
}
#'
#'@export
generateTimepoints.rSIR <- function(samples){
  minUseful <- min(which(samples[[round(length(samples)/2)]]@R != 0))
  lengthUseful <- length(samples[[round(length(samples)/2)]]@R) - minUseful
  timePoints <- c(sample(minUseful:round(minUseful + lengthUseful/3), 1),
                  sample(round(minUseful + lengthUseful/3):round(minUseful + lengthUseful*2/3), 1),
                  sample(round(minUseful + lengthUseful*2/3):length(samples[[round(length(samples)/2)]]@newR), 1))
  return(timePoints)
}
