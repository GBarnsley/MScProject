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
  RWM <- nimbleFunction(
    contains = sampler_BASE,
    setup = sampler_RW_block_TRACER_setup,
    run = sampler_RW_block_TRACER_run,
    methods = sampler_RW_block_TRACER_methods
  )
  UseMethod("buildMCMCInternal")
}
#'
#'@export
buildMCMCInternal.SIR <- function(epiModel, hyperParameters){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  control <- buildControl(epiModel, hyperParameters, "RWM" , 1)
  output$addSampler(target = c('Beta', 'Gamma'), type = RWM, control = control)
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
  control <- buildControl(epiModel, hyperParameters, "RWM" , 1)
  output$addSampler(target = c('Beta', 'Gamma'), type = RWM, control = control)
  control <- buildControl(epiModel, hyperParameters, "N+-Delta", 2)
  output$addSampler(target = "newI",
                    type = sampler,
                    control = control)
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
  control <- buildControl(epiModel, hyperParameters, "RWM" , 1)
  output$addSampler(target = c('Beta', 'Gamma'), type = RWM, control = control)
  control <- buildControl(epiModel, hyperParameters, "N+-Delta", 2)
  output$addSampler(target = "newR",
                    type = sampler,
                    control = control)
  output$addMonitors(c('Beta', 'Gamma', 'newR'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE)
  return(output)
}
#'
#'@export
buildControl <- function(epiModel, hyperParameters, algorithm, column){
  control <- hyperParameters[[algorithm]]
  control$Column <- column
  return(control)
}
