#'
#'@export
generateTimepoints <- function(epiModel){
  UseMethod("generateTimepoints")
}
#'
#'@export
generateTimepoints.SIR <- function(epiModel){
  usefulNewILength <- 0
  for(i in 1:nrow(epiModel@Samples)){
    newILength <- max(which(test@Samples[i,-c(1,2)] != 0))
    if(newILength > usefulNewILength){
      usefulNewILength <- newILength
    }
  }
  timePoints <- c(sample(1:round(usefulNewILength/3), 1),
                  sample(round(usefulNewILength/3):round(usefulNewILength*2/3), 1),
                  sample(round(usefulNewILength*2/3):usefulNewILength, 1))
  return(timePoints)
}
