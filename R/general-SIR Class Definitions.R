SIRclass <- setClass(
  "SIR",
  slots = c(
    Model = "ANY",
    MCMC = "ANY",
    Samples = "ANY",
    Metrics = "ANY"
    )
)
newISIRclass <- setClass(
  "iSIR",
  contains = "SIR"
)
newRSIRclass <- setClass(
  "rSIR",
  contains = "SIR"
)
#' Creates an object of class SIR, iSIR or rSIR
#'
#' Takes any combination of the values in the SIR model and uses dataGen to generate
#' unspecified data where it can. If S, I and newI but R and newR are then it will
#' generate an object of class iSIR that inherits from the SIR class. This class functions
#' the same for most functions other than in the metropolis hasting algorithm where it treats
#' newI as a parameter to generate. if R, I and newR are missing then it will make an object of
#' class rSIR which does the same but treats newR as a parameter.
#'
#' @param S Time series of the Susceptible Population
#' @param I Time series of the Infectious Population
#' @param R Time series of the Recovered Population
#' @param newI Time series of new infections
#' @param newR Time series of new recoveries
#' @param N The size of the population
#' @param Beta The infectious weight of one person in I
#' @param Gamma The recovery rate
#' @param t.step The amount of time in each step of the model, actual value depends on the units of Beta and Gamma
#' @return Object of class SIR
#' @export
SIR <- function(S = NULL,
                I = NULL,
                R = NULL,
                newI = NULL,
                newR = NULL,
                N = NULL,
                Beta = NULL,
                Gamma = NULL,
                t.step = 1){
  #calculating initial values from given dataS
  if(is.null(N)){
    print("Error: N must be specified")
    return(NA)
  }
  if(is.null(S)&!is.null(newI)){
    S <- rep(N-1, length(newI) + 1)
    for(i in 2:length(S)){
      S[i] <- N - sum(newI[1:i-1])
    }
  }
  if(is.null(R)&!is.null(newR)){
    R <- rep(0, length(newR) + 1)
    for(i in 2:(length(R)-1)){
      R[i] <- sum(newR[1:i-1])
    }
  }
  if(is.null(I)&!is.null(S)&!is.null(R)){
    I <- N - S - R
  }
  if(is.null(newI)&!is.null(S)){
    newI <- -diff(S)
  }
  if(is.null(newR)&!is.null(R)){
    newR <- diff(R)
  }
  #Setting up models based on nulls
  if(is.null(newI)){
    tempCode <- nimble::nimbleCode({
      # Set priors
      Beta ~ dgamma(shape = BetaShape, rate = BetaRate)
      Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
      # likelihood
      S[1] <- Pop - 1
      I[1] <- 1
      for(i in 1:TimePeriod){
        newI[i] ~ dbinom(size = S[i],
                         prob =  probGen(I[i]*Beta*t.step))
        newR[i] ~ dbinom(size = I[i], prob =  probGen(Gamma*t.step))
        S[i+1] <- S[i] - newI[i]
        I[i+1] <- I[i] + newI[i] - newR[i]
      }
    })
    return(newISIRclass(
      Model = nimble::compileNimble(
        nimble::nimbleModel(
          code = tempCode,
          constants = list(TimePeriod = length(newR)),
          data = list(newR = newR,
                      t.step = t.step,
                      Pop = N,
                      BetaShape = 1,
                      BetaRate = 1,
                      GammaShape = 1,
                      GammaRate = 1),
          inits = list(Beta = 1,
                       Gamma = 1,
                       newI = rep(0, length(newR))),
          calculate = FALSE
        )
      ),
      MCMC = NA,
      Samples = NA,
      Metrics = NA
    )
    )
  }
  else if(is.null(newR)){
    tempCode <- nimble::nimbleCode({
      # Set priors
      Beta ~ dgamma(shape = BetaShape, rate = BetaRate)
      Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
      # likelihood
      S[1] <- Pop - 1
      I[1] <- 1
      for(i in 1:TimePeriod){
        newI[i] ~ dbinom(size = S[i],
                         prob =  probGen(I[i]*Beta*t.step))
        newR[i] ~ dbinom(size = I[i], prob =  probGen(Gamma*t.step))
        S[i+1] <- S[i] - newI[i]
        I[i+1] <- I[i] + newI[i] - newR[i]
      }
    })
    return(newRSIRclass(
      Model = nimble::compileNimble(
        nimble::nimbleModel(
          code = tempCode,
          constants = list(TimePeriod = length(newI)),
          data = list(newI = newI,
                      Pop = N,
                      t.step = t.step,
                      BetaShape = 1,
                      BetaRate = 1,
                      GammaShape = 1,
                      GammaRate = 1),
          inits = list(Beta = 1,
                       Gamma = 1,
                       newR = rep(0, length(newI))),
          calculate = FALSE
        )
      ),
      MCMC = NA,
      Samples = NA,
      Metrics = NA
    )
    )
  }
  else{
    tempCode <- nimble::nimbleCode({
      # Set priors
      Beta ~ dgamma(shape = BetaShape, rate = BetaRate)
      Gamma ~ dgamma(shape = GammaShape, rate = GammaRate)
      # likelihood
      for(i in 1:TimePeriod){
        newI[i] ~ dbinom(size = S[i],
                         prob =  probGen(I[i]*Beta*t.step))
        newR[i] ~ dbinom(size = I[i], prob =  probGen(Gamma*t.step))
      }
    })
    return(SIRclass(
      Model = nimble::compileNimble(
        nimble::nimbleModel(
          code = tempCode,
          constants = list(TimePeriod = length(newI)),
          data = list(I = I,
                      S = S,
                      newI = newI,
                      newR = newR,
                      t.step = t.step,
                      BetaShape = 1,
                      BetaRate = 1,
                      GammaShape = 1,
                      GammaRate = 1),
          inits = list(Beta = 1,
                       Gamma = 1),
          calculate = FALSE
        )
      ),
      MCMC = NA,
      Samples = NA,
      Metrics = NA
    )
    )
  }
}
