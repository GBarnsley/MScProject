#' function that sets up initial values for MH algorithm.
#' @export
initialValues <- function(epiModel, hyperParameters){
  UseMethod("initialValues", epiModel)
}
#' Intial values for SIR model. Just sets Beta and Gamma to
#' the start value specified in the hyperparameters.
#' from it.
#' @export
initialValues.SIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Alpha
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Beta
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Alpha
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Beta
  return(
    epiModel
  )
}
#' Intial values for SIR model. Just sets Beta and Gamma to
#' the start value specified in the hyperparameters and then
#' generates newI.
#' from it.
#' @export
initialValues.iSIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Alpha
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Beta
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Alpha
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Beta
  timePeriod <- length(epiModel@Model$newR)
  epidemicSize <- sum(epiModel@Model$newR)
  prob2 <- probGen(
    epiModel@Model$Gamma * epiModel@Model$t.step
  )
  epiModel@Model$S[1] <- epiModel@Model$Pop - 1
  epiModel@Model$I[1] <- 1
  for(i in 1:(timePeriod-1)){
    if(epiModel@Model$I[i] != 0 & epiModel@Model$S[i] != 0){
      values <- max(epiModel@Model$newR[i] - epiModel@Model$I[i] + 1, epiModel@Model$newR[i+1] + epiModel@Model$newR[i] - epiModel@Model$I[i] + 1, 0):(epiModel@Model$S[i] - (epiModel@Model$Pop - epidemicSize))
      if(length(values) > 1){
        densities <- rep(NA, length(values))
        prob <-  probGen(
          epiModel@Model$I[i]*epiModel@Model$Beta*epiModel@Model$t.step/(epiModel@Model$Pop^epiModel@Model$Frequency)
        )
        incompleteSize <- epiModel@Model$I[i] - epiModel@Model$newR[i]
        for(j in 1:length(values)){
          densities[j] <- dbinom(values[j], epiModel@Model$S[i], prob, log = TRUE) +
            dbinom(epiModel@Model$newR[i+1], incompleteSize + values[j], prob2, log = TRUE)
        }
        prob <- densityToProbability(densities)
        epiModel@Model$newI[i] <- sample(values, 1, FALSE, prob = prob)
      }
      else{
        epiModel@Model$newI[i] <- values[1]
      }
    }
    else{
      epiModel@Model$newI[i] <- 0
    }
    epiModel@Model$S[i+1] <- epiModel@Model$S[i] - epiModel@Model$newI[i]
    epiModel@Model$I[i+1] <- epiModel@Model$I[i] + epiModel@Model$newI[i] - epiModel@Model$newR[i]
  }
  epiModel@Model$newI[timePeriod] <- 0
  epiModel@Model$S[timePeriod + 1] <- epiModel@Model$S[timePeriod]
  epiModel@Model$I[timePeriod + 1] <- epiModel@Model$I[timePeriod] - epiModel@Model$newR[i]
  return(
    epiModel
  )
}
#' Intial values for SIR model. Just sets Beta and Gamma to
#' the start value specified in the hyperparameters and then
#' generates newR.
#' from it.
#' @export
initialValues.rSIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Alpha
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Beta
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Alpha
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Beta
  timePeriod <- length(epiModel@Model$newI)
  epidemicSize <- sum(epiModel$newI) + 1
  epiModel@Model$I[1] <- 1
  for(i in 1:timePeriod){
    epiModel@Model$newR[i] <- rbinom(n = 1,
                                            size = epiModel@Model$I[i],
                                            prob = probGen(
                                              epiModel@Model$Gamma*epiModel@Model$t.step
                                            )
    )
    epiModel@Model$I[i+1] <- epiModel@Model$I[i] + epiModel@Model$newI[i] - epiModel@Model$newR[i]
  }
  return(
    epiModel
  )
}
#'
#'@export
densityToProbability <- function(densities, UB = 700, LB = 700){
  maxValue <- max(densities)
  scale <- UB - maxValue
  densities <- densities + scale
  prob <- rep(0, length(densities))
  prob[densities >= LB] <- exp(densities[densities >= LB])
  return(prob)
}
