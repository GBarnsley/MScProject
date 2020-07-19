#' Generates Many Simulated Epidemics for the SIR class
#' @param t.max The maximum length of time to run for, if NA then it will run until I=0
#' @export
simulateSIRs <- function(Beta, Gamma, Pop, N, t.step = 1, t.max = NA, seed = NA, Density = FALSE){
  if(!is.na(seed)){
    set.seed(seed)
  }
  if(identical(t.max, NA)){
    t.max <- N*10 #just a large value, there's probably a better way of doing this
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
  for(i in 1:N){
    tempEpidemic$simulate()
    maxLength <- max(max(which(tempEpidemic$newR!=0)), max(which(tempEpidemic$newI!=0), 0))
    output$newR[[i]] <- tempEpidemic$newR
    output$newI[[i]] <- tempEpidemic$newI
  }
  return(output)
}
