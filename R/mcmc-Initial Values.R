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
  for(i in 1:(length(epiModel@Model$newR)-1)){
    epiModel@Model$newI[i] <- epiModel@Model$newR[i+1]
  }
  if(epiModel@Model$newR[1] != 0){
    epiModel@Model$newI[1] <- epiModel@Model$newI[1] + epiModel@Model$newR[1]
  }
  if(sum(epiModel@Model$newI) == epiModel@Model$Pop){
    epiModel@Model$newI[max(which(epiModel@Model$newI!=0))] <-
      epiModel@Model$newI[max(which(epiModel@Model$newI!=0))] - 1
  }
  epiModel@Model$calculate()
  epiModel@Model$newI[length(epiModel@Model$newR)] <- rbinom(1,
                                                             epiModel@Model$S[length(epiModel@Model$newR)],
                                                             probGen(epiModel@Model$t.step*
                                                                       epiModel@Model$I[length(epiModel@Model$newR)]*
                                                                       epiModel@Model$Beta)
                                                             ) #THis part needs some work!! probably add parameter in to be size of epidemic etc. + prior?
  epiModel@Model$calculate()
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
  epiModel@Model$newR[1] <- 0
  for(i in 2:epiModel@Model$TimePeriod){
    epiModel@Model$newR[i] <- epiModel@Model$newI[i-1]
  }
  epiModel@Model$calculate()
  return(
    epiModel
  )
}
