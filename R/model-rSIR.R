#' Method to initialize the rSIR model, sets Beta/Gamma by their given initial
#' values then runs the NpmDelta algorithm on newR.
#' @param epiModel An object of the rSIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return rSIR class object with the initial values
#' @export
initialValues.rSIR <- function(epiModel, hyperParameters){
  epiModel@Model$Beta <- hyperParameters$`Initial Values`$Beta
  epiModel@Model$Gamma <- hyperParameters$`Initial Values`$Gamma
  epiModel@Model$BetaShape <- hyperParameters$Priors$Beta$Shape
  epiModel@Model$BetaRate <- hyperParameters$Priors$Beta$Rate
  epiModel@Model$GammaShape <- hyperParameters$Priors$Gamma$Shape
  epiModel@Model$GammaRate <- hyperParameters$Priors$Gamma$Rate
  epiModel@Model$newR <- rep(0, length(epiModel@Model$newI))
  epiModel@Model$newR[length(epiModel@Model$newI)] <- sum(epiModel@Model$newI) + 1
  mcmc <- configureMCMC(epiModel@Model, nodes = NULL)
  mcmc$addSampler(target = "newR",
                  type = sampler,
                  control = list(
                    TMax = 20,
                    DeltaMax = 20,
                    R = hyperParameters$`Initial Values`$Runs
                  ))
  mcmc <- buildMCMC(
    mcmc
  )
  mcmc <- compileNimble(mcmc, project = epiModel@Model, resetFunctions = TRUE)
  mcmc$run(1)
  return(
    epiModel
  )
}
#' Builds the MCMC for the rSIR model. Applies a block-RWM to Beta and Gamma and
#' the NpmDelta sampler to newR.
#' @param epiModel An object of the rSIR class
#' @param hyperParameters A list of lists of the hyper-parameters for the epidemic model and MCMC
#' @return rSIR class with a compiled MCMC
#' @export
buildMCMCInternal.rSIR <- function(epiModel, hyperParameters){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = "newR",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addMonitors(c('Beta', 'Gamma', 'newR'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE)
  return(output)
}
