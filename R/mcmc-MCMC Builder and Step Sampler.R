#'
#'@export
buildMCMCInternal <- function(epiModel, hyperParameters){
  if(identical(hyperParameters$`True Values`, NA)|identical(hyperParameters$`True Values`, NULL)){
    stepSampler <- nimbleFunction(
      contains = sampler_BASE,
      setup = function(model, mvSaved, target, control) {
        minStep <- as.integer(control$`Minimum Step Size`)
        maxStep <- as.integer(control$`Maximum Step Size`)
        maxChange <- as.integer(control$`Maximum Change`)
        calcNodes <- model$getDependencies(target)
        evalCol <- as.integer(control$`Assigned Column`)
      },
      ## the run function
      run = function() {
        positions <- which(model[[target]]!=0)
        pos <- rcat(n = 1, prob = rep(1/length(positions), length(positions)))
        position <- positions[pos]

        directions <- c(-1,1)
        pos <- rcat(n = 1, prob = rep(1/length(directions), length(directions)))
        direction <- directions[pos]

        sizes <- minStep:maxStep
        pos <- rcat(n = 1, prob = rep(1/length(sizes), length(sizes)))
        size <- sizes[pos]

        amounts <- min(maxChange, model[[target]][position])
        amount <- rcat(n = 1, prob = rep(1/amounts, amounts))

        newPosition <- position + direction*size

        if(newPosition > 0 & newPosition <= length(model[[target]])){
          model_lp_initial <- getLogProb(model, calcNodes)

          logProbForward <- -log(sum(model[[target]]!=0)) -
            log(min(maxChange, model[[target]][position]))
          model[[target]][position] <<- model[[target]][position] - amount
          model[[target]][newPosition] <<- model[[target]][newPosition] + amount
          logProbBackward <- -log(sum(model[[target]]!=0)) -
            log(min(maxChange, model[[target]][newPosition]))

          model_lp_proposed <- calculate(model, calcNodes)

          log_MH_ratio <- model_lp_proposed - model_lp_initial + logProbBackward - logProbForward

          u <- runif(1, 0, 1)
          if(u < exp(log_MH_ratio)){
            model[["tracers"]][1:3,evalCol] <<- c(model[["tracers"]][1,evalCol] + 1,
                                                  model[["tracers"]][2,evalCol] + size^2,
                                                  0)
            jump <- TRUE
          }
          else{
            jump <- FALSE
          }
        }
        else{
          jump <- FALSE
        }
        ## if we accepted the proposal, then store the updated
        ## values and logProbs from 'model' into 'mvSaved'.
        ## if the proposal was not accepted, restore the values
        ## and logProbs from 'mvSaved' back into 'model'.
        if(jump) copy(from = model, to = mvSaved, row = 1,
                      nodes = calcNodes, logProb = TRUE)
        else copy(from = mvSaved, to = model, row = 1,
                  nodes = calcNodes, logProb = TRUE)
      },
      methods = list(
        reset = function () {}
      )
    )
  }
  else{
    stepSampler <- nimbleFunction(
      contains = sampler_BASE,
      setup = function(model, mvSaved, target, control) {
        minStep <- as.integer(control$`Minimum Step Size`)
        maxStep <- as.integer(control$`Maximum Step Size`)
        maxChange <- as.integer(control$`Maximum Change`)
        calcNodes <- model$getDependencies(target)
        evalCol <- as.integer(control$`Assigned Column`)
        trueValue <- control$`True Value`
      },
      ## the run function
      run = function() {
        positions <- which(model[[target]]!=0)
        pos <- rcat(n = 1, prob = rep(1/length(positions), length(positions)))
        position <- positions[pos]

        directions <- c(-1,1)
        pos <- rcat(n = 1, prob = rep(1/length(directions), length(directions)))
        direction <- directions[pos]

        sizes <- minStep:maxStep
        pos <- rcat(n = 1, prob = rep(1/length(sizes), length(sizes)))
        size <- sizes[pos]

        amounts <- min(maxChange, model[[target]][position])
        amount <- rcat(n = 1, prob = rep(1/amounts, amounts))

        newPosition <- position + direction*size

        if(newPosition > 0 & newPosition <= length(model[[target]])){
          model_lp_initial <- getLogProb(model, calcNodes)

          logProbForward <- -log(sum(model[[target]]!=0)) -
            log(min(maxChange, model[[target]][position]))
          model[[target]][position] <<- model[[target]][position] - amount
          model[[target]][newPosition] <<- model[[target]][newPosition] + amount
          logProbBackward <- -log(sum(model[[target]]!=0)) -
            log(min(maxChange, model[[target]][newPosition]))

          model_lp_proposed <- calculate(model, calcNodes)

          log_MH_ratio <- model_lp_proposed - model_lp_initial + logProbBackward - logProbForward

          u <- runif(1, 0, 1)
          if(u < exp(log_MH_ratio)){
            model[["tracers"]][1:3,evalCol] <<- c(model[["tracers"]][1,evalCol] + 1,
                                                  model[["tracers"]][2,evalCol] + size^2,
                                                  model[["tracers"]][3,evalCol] +
                                                  mean((trueValue - model[[target]])^2))
            jump <- TRUE
          }
          else{
            jump <- FALSE
          }
        }
        else{
          jump <- FALSE
        }
        ## if we accepted the proposal, then store the updated
        ## values and logProbs from 'model' into 'mvSaved'.
        ## if the proposal was not accepted, restore the values
        ## and logProbs from 'mvSaved' back into 'model'.
        if(jump) copy(from = model, to = mvSaved, row = 1,
                      nodes = calcNodes, logProb = TRUE)
        else copy(from = mvSaved, to = model, row = 1,
                  nodes = calcNodes, logProb = TRUE)
      },
      methods = list(
        reset = function () {}
      )
    )
  }
  UseMethod("buildMCMCInternal")
}
#'
#'@export
buildMCMCInternal.SIR <- function(epiModel, hyperParameters){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = 'Beta', type = 'RW', control = list(scale = hyperParameters$`Random Walk`$Beta, adaptive = FALSE))
  output$addSampler(target = 'Gamma', type = 'RW', control = list(scale = hyperParameters$`Random Walk`$Gamma, adaptive = FALSE))
  output$addMonitors(c('Beta', 'Gamma'))
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
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = 'Beta', type = 'RW', control = list(scale = hyperParameters$`Random Walk`$Beta, adaptive = FALSE))
  output$addSampler(target = 'Gamma', type = 'RW', control = list(scale = hyperParameters$`Random Walk`$Gamma, adaptive = FALSE))
  control <- hyperParameters$`Step Proposal`
  control$`Assigned Column` <- 3
  control$`True Value` <- hyperParameters$`True Values`$newI
  output$addSampler(target = "newI",
                    type = stepSampler,
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
  output$addSampler(target = 'Beta', type = 'RW', control = list(scale = hyperParameters$`Random Walk`$Beta, adaptive = FALSE))
  output$addSampler(target = 'Gamma', type = 'RW', control = list(scale = hyperParameters$`Random Walk`$Gamma, adaptive = FALSE))
  control <- hyperParameters$`Step Proposal`
  control$`Assigned Column` <- 3
  control$`True Value` <- hyperParameters$`True Values`$newR
  output$addSampler(target = "newR",
                    type = stepSampler,
                    control = control)
  output$addMonitors(c('Beta', 'Gamma', 'newIR'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE)
  return(output)
}
