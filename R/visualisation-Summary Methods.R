#'A generic function to return a summary model using either the mean, median or mode of the estimated paratmers
#'Unused at the moment.
#'@export
summarise <- function(samples, average = mean){
 UseMethod("summarise", samples)
}
#' Method that returns a epidemic model with mean values of its
#' parameters, calculated from the samples, for the SIR class
#' @export
summarise.SIR <- function(samples, average){
  output <- samples[[1]]
  Betas <- rep(NA, length(samples))
  Gammas <- rep(NA,length(samples))
  for(i in 1:length(Betas)){
    Betas[i] <- samples[[i]]@Beta
    Gammas[i] <- samples[[i]]@Gamma
  }
  output@Beta <- average(Betas)
  output@Gamma <- average(Gammas)
  return(output)
}
#' Method that returns a epidemic model with average values of its
#' parameters, calculated from the samples, for the iSIR class
#' @export
summarise.iSIR <- function(samples, average){
  output <- samples[[1]]
  Betas <- rep(NA, length(samples))
  Gammas <- rep(NA, length(samples))
  newIs <- matrix(NA,
                 nrow = length(samples),
                 ncol = length(samples[[1]]@newI))
  Ss <- matrix(NA,
              nrow = length(samples),
              ncol = length(samples[[1]]@S))
  Is <- matrix(NA,
              nrow = length(samples),
              ncol = length(samples[[1]]@I))
  for(i in 1:length(Betas)){
    Betas[i] <- samples[[i]]@Beta
    Gammas[i] <- samples[[i]]@Gamma
    for(j in 1:length(samples[[1]]@S)){
      if(j < length(samples[[1]]@S)){
        newIs[i,j] <- samples[[i]]@newI[j]
      }
      Ss[i,j] <- samples[[i]]@S[j]
      Is[i,j] <- samples[[i]]@I[j]
    }
  }
  S <- rep(NA, length(samples[[1]]@S))
  I <- rep(NA, length(samples[[1]]@I))
  newI <- rep(NA, length(samples[[1]]@newI))
  for(i in 1:length(samples[[1]]@S)){
    S[i] <- average(Ss[,i])
    I[i] <- average(Is[,i])
    if(i < length(samples[[1]]@S)){
      newI <- average(newIs[,i])
    }
  }
  output@Beta <- average(Betas)
  output@Gamma <- average(Gammas)
  output@S <- S
  output@I <- I
  output@newI <- newI
  return(output)
}
#' Method that returns a epidemic model with mean values of its
#' parameters, calculated from the samples, for the rSIR class
#' @export
summarise.rSIR <- function(samples, average){
  output <- samples[[1]]
  Betas <- rep(NA, length(samples))
  Gammas <- rep(NA, length(samples))
  newRs <- matrix(NA,
                  nrow = length(samples),
                  ncol = length(samples[[1]]@newR))
  Rs <- matrix(NA,
               nrow = length(samples),
               ncol = length(samples[[1]]@R))
  Is <- matrix(NA,
               nrow = length(samples),
               ncol = length(samples[[1]]@I))
  for(i in 1:length(Betas)){
    Betas[i] <- samples[[i]]@Beta
    Gammas[i] <- samples[[i]]@Gamma
    for(j in 1:length(samples[[1]]@R)){
      if(j < length(samples[[1]]@R)){
        newRs[i,j] <- samples[[i]]@newR[j]
      }
      Rs[i,j] <- samples[[i]]@R[j]
      Is[i,j] <- samples[[i]]@I[j]
    }
  }
  R <- rep(NA, length(samples[[1]]@R))
  I <- rep(NA, length(samples[[1]]@I))
  newR <- rep(NA, length(samples[[1]]@newR))
  for(i in 1:length(samples[[1]]@R)){
    R[i] <- average(Rs[,i])
    I[i] <- average(Is[,i])
    if(i < length(samples[[1]]@S)){
      newR <- average(newRs[,i])
    }
  }
  output@Beta <- average(Betas)
  output@Gamma <- average(Gammas)
  output@R <- R
  output@I <- I
  output@newR <- newR
  return(output)
}
