#' Generates Many Simulated Epidemics for the SIR class
#' @export
simulateSIRs <- function(Beta, Gamma, Pop, N, t.step = 1, t.max = 100, seed = NA){
  output <- list(training = list(),
                 test = list())
  if(!is.na(seed)){
    set.seed(seed)
  }
  tempCode <- nimble::nimbleCode({
    # likelihood
    S[1] <- Pop - 1
    I[1] <- 1
    for(i in 1:TimePeriod){
      newI[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Beta*t.step))
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
                  t.step = t.step)
    )
  )
  for(i in 1:N){
    tempEpidemic$simulate()
    output$newR[[i]] <- tempEpidemic$newR
    output$newI[[i]] <- tempEpidemic$newI
  }
  return(output)
}
