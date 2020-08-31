#' Function to simulate SIR epidemics.
#' Automatically discounts failed epidemics (if the epidemic fails to infect more
#' than a 100th of the total population) and epidemics that fail to finish before
#' the maximum time period. The number of failures is printed at the end.
#' @param Beta The infective weight of a single infectious individual of the epidemic
#' @param Gamma The recovery rate of the epidemic
#' @param Pop Total population in the epidemic
#' @param N The number of epidemics to simulate
#' @param t.step The time-step the model follows, only affects the scale of Beta and Gamma
#' @param seed The random seed for the start of the simulations
#' @param Frequency A true/false value, specifying if the model has frequency-based transmission
#' @param t.max The maximum length of time to run for, Ideally should be adjusted so that the number of incomplete epidemics is close to 0
#' @return A list of transition vectors for infections and recoveries
#' @export
simulateSIRs <- function(Beta, Gamma, Pop, N, t.step = 1, t.max = NA, seed = NA, Frequency = TRUE){
  if(!is.na(seed)){
    set.seed(seed)
  }
    tempCode <- nimble::nimbleCode({
      # likelihood
      S[1] <- Pop - 1
      I[1] <- 1
      for(i in 1:TimePeriod){
        newI[i] ~ dbinom(size = S[i],
                         prob =  probGen(I[i]*Beta*t.step/(Pop^Frequency)))
        newR[i] ~ dbinom(size = I[i],
                         prob =  probGen(Gamma*t.step))
        S[i+1] <- S[i] - newI[i]
        I[i+1] <- I[i] + newI[i] - newR[i]
      }
    })
  tempEpidemic <- nimble::compileNimble(
    nimble::nimbleModel(
      code = tempCode,
      constants = list(TimePeriod = t.max),
      data = list(Pop = Pop,
                  Beta = Beta,
                  Gamma = Gamma,
                  t.step = t.step,
                  Frequency = Frequency*1)
    ))
  output <- list(newI = vector("list", N),
                 newR = vector("list", N))
  unfinished <- 0
  failed <- 0
  successfull <- 0
  while(successfull < N){
    if(failed + unfinished >= 10*N){
      print("Failed to generate values in resonable time, please adjust parameters")
      break
    }
    tempEpidemic$simulate()
    test <- tempEpidemic$I==0
    if(sum(tempEpidemic$newI) < Pop/100){
      failed <- failed + 1
    }
    else if(any(test)&sum(test)<t.max-1){
      successfull <- successfull + 1
      end <- min(which(tempEpidemic$I==0)) - 1
      output$newR[[successfull]] <- tempEpidemic$newR[1:end]
      output$newI[[successfull]] <- tempEpidemic$newI[1:end]
    }
    else{
      unfinished <- unfinished + 1
    }
  }
  print("Unfinished epidemics:")
  print(unfinished)
  print("Failed epidemics:")
  print(failed)
  return(output)
}
