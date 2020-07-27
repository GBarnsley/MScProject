#' Generates Many Simulated Epidemics for the SIR class
#' @param t.max The maximum length of time to run for, Ideally should be adjusted so that the number of incomplete epidemics is close to 0
#' @export
simulateSIRs <- function(Beta, Gamma, Pop, N, t.step = 1, t.max = NA, seed = NA, Density = FALSE){
  if(!is.na(seed)){
    set.seed(seed)
  }
    tempCode <- nimble::nimbleCode({
      # likelihood
      S[1] <- Pop - 1
      I[1] <- 1
      for(i in 1:TimePeriod){
        newI[i] ~ dbinom(size = S[i],
                         prob =  probGen(I[i]*Beta*t.step/(Pop^Density)))
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
                  Density = Density*1)
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
