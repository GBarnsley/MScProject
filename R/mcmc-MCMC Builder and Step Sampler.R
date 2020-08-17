#'
#'@export
buildMCMCInternal <- function(epiModel, hyperParameters){
  sampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = stepSampler_setup,
    run = stepSampler_run,
    methods = list(
      reset = function() {}
    )
  )
  UseMethod("buildMCMCInternal")
}
#'
#'@export
buildMCMCInternal.SIR <- function(epiModel, hyperParameters){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  return(compileNimble(
    buildMCMC(
      output
    ),
    project = epiModel@Model,
    resetFunctions = TRUE
  ))
}
#'
#'@export
buildMCMCInternal.iSIR <- function(epiModel, hyperParameters){
  if(!is.null(hyperParameters$newI$Bounded)){
    if(hyperParameters$newI$Bounded){
      sampler <- nimbleFunction(
        contains = sampler_BASE,
        setup = stepSampler_setup,
        run = stepSampler_run_bounded,
        methods = list(
          reset = function() {}
        )
      )
    }}
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Beta', 'Gamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addMonitors(c('Beta', 'Gamma', 'newI'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE)
  return(output)
}
#'
#'@export
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

