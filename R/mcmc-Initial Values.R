#' function that sets up initial values for MH algorithm.
#' @export
initialValues <- function(epiModel, hyperParameters){
  sampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = stepSampler_setup,
    run = stepSampler_run,
    methods = list(
      reset = function() {}
    )
  )
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
  epiModel@Model$newI <- rep(0, length(epiModel@Model$newR))
  epiModel@Model$newI[1] <- sum(epiModel@Model$newR) - 1
  mcmc <- configureMCMC(epiModel@Model, nodes = NULL)
  mcmc$addSampler(target = "newI",
                    type = sampler,
                    control = list(
                      TMax = 20,
                      DeltaMax = 20,
                      R = hyperParameters$`Initial Values`$Runs,
                      Column = 1
                    ))
  mcmc <- buildMCMC(
    mcmc
  )
  mcmc <- compileNimble(mcmc, project = epiModel@Model, resetFunctions = TRUE)
  mcmc$run(1)
  epiModel@Model$tracers <- matrix(0, nrow = 2, ncol = 2)
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
  epiModel@Model$newR <- rep(0, length(epiModel@Model$newI))
  epiModel@Model$newR[length(epiModel@Model$newI)] <- sum(epiModel@Model$newI) + 1
  mcmc <- configureMCMC(epiModel@Model, nodes = NULL)
  mcmc$addSampler(target = "newR",
                  type = sampler,
                  control = list(
                    TMax = 20,
                    DeltaMax = 20,
                    R = hyperParameters$`Initial Values`$Runs,
                    Column = 1
                  ))
  mcmc <- buildMCMC(
    mcmc
  )
  mcmc <- compileNimble(mcmc, project = epiModel@Model, resetFunctions = TRUE)
  mcmc$run(1)
  epiModel@Model$tracers <- matrix(0, nrow = 2, ncol = 2)
  return(
    epiModel
  )
}
