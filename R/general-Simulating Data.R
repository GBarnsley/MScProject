#' Generates Simulated Epidemics
#'
#' Given an model with specified parameters, e.g. Beta, Gamma and population size,
#' this function returns a stochastically simulated epidemic. The model used will
#' be of the same type as the provided model. Currently only accepts SIR models.
#'
#' @param epiModel An epidemic model with know parameters
#' @param t.max the largest number of time-steps to be taken
#' @return an object containing the generated data
#' @export
simulate.SIR <- function(epiModel, t.max = 1000){
  if(is.na(epiModel@Beta)|is.na(epiModel@Gamma)|is.na(epiModel@N)|is.na(epiModel@t.step)){
    print("Missing parameter, simulation cannot complete")
  }
  else{
    S <- epiModel@N-1; I <- 1; R <- 0 #initial conditions
    T <- seq(0, t.max, by = epiModel@t.step) #sequence of steps
    res <- matrix(nrow = length(T), ncol = 6) #matrix to hold results
    colnames(res) = c('time','S','I','R', 'newI', 'newR')
    res[1,] = c(0, S, I, R, NA, NA) #first step, initial conditions
    for (step in 2:length(T)){
      #calculating new infections
      Pinf <- probGen(epiModel@Beta*I*epiModel@t.step) #pinf for time t
      Inew <- rbinom(1, size = S, prob = Pinf) #draw a sample from binomial using Pinf
      #calculating removals
      Prem <- probGen(epiModel@Gamma*epiModel@t.step) #prem for time t
      Rnew <- rbinom(1, size = I, prob = Prem) #draw another binomial sample
      #updating states
      S <- S - Inew
      I <- I + Inew - Rnew
      R <- R + Rnew
      #storing states
      res[step,] <- c(T[step], S, I, R, NA, NA)
      res[step-1,5] <- Inew
      res[step-1,6] <- Rnew
      if(I == 0){
        break
      }
    }
    return(SIR(S = na.omit(res[,2]),
               I = na.omit(res[,3]),
               R = na.omit(res[,4]),
               newI = na.omit(res[,5]),
               newR = na.omit(res[,6]),
               Beta = epiModel@Beta,
               Gamma = epiModel@Gamma,
               N = epiModel@N,
               t.step = epiModel@t.step))
  }
}
#' Generates Simulated Epidemics for iSIR
#'
#' Given an model with specified parameters, e.g. Beta, Gamma and population size,
#' R and newR this function returns a stochastically simulated epidemic.
#'
#' @param epiModel An epidemic model with know parameters
#' @return an object containing the generated data
#' @export
simulate.iSIR <- function(epiModel){
  if(epiModel@Beta <= 0|epiModel@Gamma <= 0){
    return(NA)
  }
  epiModel@S <- c(epiModel@N) - 1
  epiModel@I <- c(1)
  for(t in 1:length(epiModel@newR)){
    if(epiModel@S[t] > 0){
      epiModel@newI[t] <- rbinom(1, epiModel@S[t], probGen(epiModel@Beta*epiModel@t.step*epiModel@I[t]))
    }
    else{
      epiModel@newI[t] <- 0
    }
    epiModel@S[t+1] <- epiModel@S[t] - epiModel@newI[t]
    epiModel@I[t+1] <- epiModel@N - epiModel@S[t+1] - epiModel@R[t+1]
    if(epiModel@S[t+1]<0|epiModel@I[t+1]<0){
      return(NA)
    }
    #print(c(epiModel@S[t],epiModel@newI[t],epiModel@I[t], epiModel@newR[t], epiModel@R[t]))
  }
  #end <- length(epiModel@newR) + 1
  #print(c(epiModel@S[end],epiModel@newI[end],epiModel@I[end], epiModel@newR[end], epiModel@R[end]))
  return(epiModel)
}
#' Generates Simulated Epidemics for rSIR
#'
#' Given an model with specified parameters, e.g. Beta, Gamma and population size,
#' R and newR this function returns a stochastically simulated epidemic.
#'
#' @param epiModel An epidemic model with know parameters
#' @return an object containing the generated data
#' @export
simulate.rSIR <- function(epiModel){
  if(epiModel@Beta <= 0|epiModel@Gamma <= 0){
    return(NA)
  }
  epiModel@R <- c(0)
  epiModel@I <- c(1)
  for(t in 1:length(epiModel@newI)){
    epiModel@newR[t] <- rbinom(1, epiModel@I[t], probGen(epiModel@Gamma*epiModel@t.step))
    epiModel@R[t+1] <- epiModel@R[t] + epiModel@newR[t]
    epiModel@I[t+1] <- epiModel@N - epiModel@S[t+1] - epiModel@R[t+1]
    if(epiModel@I[t]<0){
      return(NA)
    }
  }
  return(epiModel)
}
#' Generates Many Simulated Epidemics
#'
#' Given an model with specified parameters, e.g. Beta, Gamma and population size,
#' this function returns 'n' stochastically simulated epidemics. Will split these into
#' testing data, with only the number of new infections, N and the timestep and into
#' the complete data. Currently only accepts SIR models.
#'
#' @param epiModel An epidemic model with know parameters
#' @param t.max the largest number of time-steps to be taken
#' @param N number of epidemics to simulate
#' @return a list containing two lists of the generated data
#' @export
simulateDataSet <- function(epiModel, N, t.max){
  UseMethod("simulateDataSet")
}
#' Generates Many Simulated Epidemics for the SIR class
#' @export
simulateDataSet.SIR <- function(epiModel, N, t.max = 1000, seed = NA){
  output <- list(test = list(),
                 complete = list())
  if(!is.na(seed)){
    set.seed(seed)
  }
  for(i in 1:N){
    tempEpidemic <- simulate(epiModel, t.max = t.max)
    output$test[[i]] <- SIR(R = tempEpidemic@R,
                          newR = tempEpidemic@newR,
                          N = tempEpidemic@N,
                          t.step = tempEpidemic@t.step)
    output$complete[[i]] <- tempEpidemic
  }
  return(output)
}
