#' Generates population data from given parameters
#'
#' Provided values of population sizes but not changes at each time step
#' this function calculates the changes and vice-versa. It will add these
#' values to a class of format provided. The purpose is that if we already
#' calculate these values here we won't have to reculate them everytime
#' we wan to find the likelihood etc. Currently on takes objects of the
#' SIR class.
#'
#' @param epiModel epidemic model in the format of a specified class
#' @return object in the same format as the given object
#' @export
dataGen <- function(epiModel){
  UseMethod(
    "dataGen"
  )
}
#' Generates population data from given parameters
#' for the SIR class.
#' @export
dataGen.SIR <- function(epiModel){
  #Generates what values it can from the supplied values
  if(any(is.null(epiModel@S))&!is.null(epiModel@newI)){
    epiModel@S <- c(epiModel@N - 1)
    for(i in 1:length(epiModel@newI)){
      epiModel@S[i+1] <- epiModel@S[i] - epiModel@newI[i]
    }
    if(any(epiModel@S < 0)){
      print("Error: negative value in data generation")
      return(NA)
    }
  }
  if(any(is.null(epiModel@R))&!is.null(epiModel@newR)){
    epiModel@R <- c()
    epiModel@R[1] <- 0
    for(i in 1:length(epiModel@newR)){
      epiModel@R[i+1] <- epiModel@R[i] + epiModel@newR[i]
    }
    if(any(epiModel@R < 0)){
      print("Error: negative value in data generation")
      return(NA)
    }
  }
  if(any(is.null(epiModel@I))&!is.null(epiModel@S)&!is.null(epiModel@R)){
    epiModel@I <- rep(epiModel@N, length(epiModel@S)) -
      epiModel@S -
      epiModel@R
    if(any(epiModel@I < 0)){
      print("Error: negative value in data generation")
      return(NA)
    }
  }
  if(any(is.null(epiModel@newI))&!is.null(epiModel@S)){
    epiModel@newI <- - diff(epiModel@S)
    if(any(epiModel@newI < 0)){
      print("Error: negative value in data generation")
      return(NA)
    }
  }
  if(any(is.null(epiModel@newR))&!is.null(epiModel@R)){
    epiModel@newR <- diff(epiModel@R)
    if(any(epiModel@newR < 0)){
      print("Error: negative value in data generation")
      return(NA)
    }
  }
  return(epiModel)
}
