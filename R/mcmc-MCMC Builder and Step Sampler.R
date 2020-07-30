#'
#'@export
buildMCMCInternal <- function(epiModel, hyperParameters){
  stepSampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      minStep <- as.integer(control$`Minimum Step Size`)
      maxStep <- as.integer(control$`Maximum Step Size`)
      maxChange <- as.integer(control$`Maximum Change`)
      calcNodes <- model$getDependencies(target)
      runs <- as.integer(control$`Runs per Random Walk`)
      evalCol <- as.integer(control$Column)
      trueValue <- control$True
    },
    ## the run function
    run = function() {
      model_lp_initial <- getLogProb(model, calcNodes)
      for(i in 1:runs){
      positions <- which(model[[target]]!=0)
      pos <- rcat(n = 1, prob = rep(1/length(positions), length(positions)))
      position <- positions[pos]

      if(position == length(model[[target]])){
        directions <- c(-1)
        direction <- -1
      }
      else if(position == 1){
        directions <- c(1)
        direction <- 1
      }
      else{
        directions <- c(-1,1)
        pos <- rcat(n = 1, prob = c(0.5,0.5))
        direction <- directions[pos]
      }

      sizes <- minStep:min(maxStep, (position - 1)*(direction == -1) +
                             (length(model[[target]]) - position)*(direction == 1))
      pos <- rcat(n = 1, prob = rep(1/length(sizes), length(sizes)))
      size <- sizes[pos]

      amounts <- min(maxChange, model[[target]][position])
      amount <- rcat(n = 1, prob = rep(1/amounts, amounts))

      newPosition <- position + direction*size

      sampler_lp_proposed <-
        (-
           #Choosing that time point
           log(sum(model[[target]]!=0)) -
           #choosing that that amount point of points to move
           log(amounts) -
           #Choosing the size of the move
           log(length(sizes)) -
           #choosing direction
           log(length(directions))
        )
      model[[target]][position] <<- model[[target]][position] - amount
      model[[target]][newPosition] <<- model[[target]][newPosition] + amount
      model_lp_proposed <- calculate(model, calcNodes)
      sampler_lp_initial <-
        (-
           #Choosing that time point
           log(sum(model[[target]]!=0)) -
           #choosing that that amount point of points to move
           log(min(maxChange, model[[target]][newPosition])) -
           #choosing the size of the move
           log(min(
             maxStep,
             (newPosition - 1)*(-direction == -1) +
               (length(model[[target]]) - newPosition)*(-direction == 1)
           ) - minStep + 1) -
           #choosing direction
           log(1 + 1*(newPosition != length(model[[target]]) & newPosition != 1))
        )


      log_MH_ratio <- (model_lp_proposed - sampler_lp_proposed) - (model_lp_initial - sampler_lp_initial)

      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        model_lp_initial <- model_lp_proposed
        model[["tracers"]][1:3,evalCol] <<- c(model[["tracers"]][1,evalCol] + 1/runs,
                                              model[["tracers"]][2,evalCol] + size^2/runs,
                                              model[["tracers"]][3,evalCol] +
                                                mean((trueValue - model[[target]])^2)/runs)
        jump <- TRUE
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
      }
    },
    methods = list(
      reset = function () {}
    )
  )
  sampler_RW_TRACER <- nimbleFunction(
    name = 'sampler_RW',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      ## control list extraction
      logScale            <- if(!is.null(control$log))                 control$log                 else FALSE
      reflective          <- if(!is.null(control$reflective))          control$reflective          else FALSE
      adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
      adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 200
      adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
      scale               <- if(!is.null(control$scale))               control$scale               else 1
      ## node list generation
      targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
      calcNodes <- model$getDependencies(target)
      calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
      ## numeric value generation
      scaleOriginal <- scale
      timesRan      <- 0
      timesAccepted <- 0
      timesAdapted  <- 0
      scaleHistory  <- c(0, 0)   ## scaleHistory
      acceptanceHistory  <- c(0, 0)   ## scaleHistory
      ## tracer values
      evalCol <- as.integer(control$Column)
      trueValue <- control$True
      if(nimbleOptions('MCMCsaveHistory')) {
        saveMCMChistory <- TRUE
      } else saveMCMChistory <- FALSE
      optimalAR     <- 0.44
      gamma1        <- 0
      ## checks
      if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
      if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
      if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
      if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
      if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
    },
    run = function() {
      currentValue <- model[[target]]
      propLogScale <- 0
      if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
      propValue <- currentValue * exp(propLogScale)
      } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
      if(reflective) {
        lower <- model$getBound(target, 'lower')
        upper <- model$getBound(target, 'upper')
        while(propValue < lower | propValue > upper) {
          if(propValue < lower) propValue <- 2*lower - propValue
          if(propValue > upper) propValue <- 2*upper - propValue
        }
      }
      model[[target]] <<- propValue
      logMHR <- calculateDiff(model, target)
      if(logMHR == -Inf) {
        nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
        ## Drawing a random number is needed during first testing
        ## of this step in order to keep the random numbers identical
        ## to old behavior to see if tests that depend on particular
        ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
        ## we call decide() to ensure same behavior.
        ## jump <- decide(logMHR)
        ## When new behavior is acceptable, we can remove the above line
        ## and uncomment the following:
        jump <- FALSE
      } else {
        logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf) + propLogScale
        jump <- decide(logMHR)
        if(jump){
          model[["tracers"]][1:3,evalCol] <<- c(model[["tracers"]][1,evalCol] + 1,
                                                model[["tracers"]][2,evalCol] +
                                                  (propValue - currentValue)^2,
                                                model[["tracers"]][3,evalCol] +
                                                  (trueValue - propValue)^2)
          nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        else{
          nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
      }
      if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
      adaptiveProcedure = function(jump = logical()) {
        timesRan <<- timesRan + 1
        if(jump)     timesAccepted <<- timesAccepted + 1
        if(timesRan %% adaptInterval == 0) {
          acceptanceRate <- timesAccepted / timesRan
          timesAdapted <<- timesAdapted + 1
          if(saveMCMChistory) {
            setSize(scaleHistory, timesAdapted)                 ## scaleHistory
            scaleHistory[timesAdapted] <<- scale                ## scaleHistory
            setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
            acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
          }
          gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
          gamma2 <- 10 * gamma1
          adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
          scale <<- scale * adaptFactor
          ## If there are upper and lower bounds, enforce a maximum scale of
          ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
          ## Otherwise, for a poorly-informed posterior,
          ## the scale could grow without bound to try to reduce
          ## acceptance probability.  This creates enormous cost of
          ## reflections.
          if(reflective) {
            lower <- model$getBound(target, 'lower')
            upper <- model$getBound(target, 'upper')
            if(scale >= 0.5*(upper-lower)) {
              scale <<- 0.5*(upper-lower)
            }
          }
          timesRan <<- 0
          timesAccepted <<- 0
        }
      },
     reset = function() {}
    )
  )
  UseMethod("buildMCMCInternal")
}
#'
#'@export
buildMCMCInternal.SIR <- function(epiModel, hyperParameters){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  control <- buildControl(epiModel, hyperParameters, "Beta", 1)
  output$addSampler(target = 'Beta', type = sampler_RW_TRACER, control = control)
  control <- buildControl(epiModel, hyperParameters, "Gamma", 2)
  output$addSampler(target = 'Gamma', type = sampler_RW_TRACER, control = control)
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
  if(!is.null(hyperParameters$newI$Bounded)){
  if(hyperParameters$newI$Bounded){
    stepSampler <- nimbleFunction(
      contains = sampler_BASE,
      setup = function(model, mvSaved, target, control) {
        minStep <- as.integer(control$`Minimum Step Size`)
        maxStep <- as.integer(control$`Maximum Step Size`)
        maxChange <- as.integer(control$`Maximum Change`)
        runs <- as.integer(control$`Runs per Random Walk`)
        calcNodes <- model$getDependencies(target)
        evalCol <- as.integer(control$Column)
        trueValue <- control$True
      },
      ## the run function
      run = function() {
        for(i in 1:runs){
        positions <- which(model[[target]]!=0)
        pos <- rcat(n = 1, prob = rep(1/length(positions), length(positions)))
        position <- positions[pos]

        if(position == length(model[[target]])){
          directions <- c(-1)
          direction <- -1
        }
        else if(position == 1){
          directions <- c(1)
          direction <- 1
        }
        else{
          directions <- c(-1,1)
          pos <- rcat(n = 1, prob = c(0.5,0.5))
          direction <- directions[pos]
        }

        sizes <- minStep:min(maxStep, (position - 1)*(direction == -1) +
                               (length(model[[target]]) - position)*(direction == 1))
        pos <- rcat(n = 1, prob = rep(1/length(sizes), length(sizes)))
        size <- sizes[pos]

        newPosition <- position + direction*size

        amounts <- min(maxChange, model[[target]][position])
        if(direction == 1){
          amounts <- min(amounts, min(model[["I"]][(position + 1):newPosition]) - 1)
        }
        amount <- rcat(n = 1, prob = rep(1/amounts, amounts))

        model_lp_initial <- getLogProb(model, calcNodes)
        model_lp_proposed <- -
          (-
             #Choosing that time point
             log(sum(model[[target]]!=0)) -
             #choosing that that amount point of points to move
             log(amounts) -
             #Choosing the size of the move
             log(length(sizes)) -
             #choosing direction
             log(length(directions))
          )
        model[[target]][position] <<- model[[target]][position] - amount
        model[[target]][newPosition] <<- model[[target]][newPosition] + amount
        model_lp_proposed <- model_lp_proposed + calculate(model, calcNodes)
        amountBackward <- min(maxChange, model[[target]][newPosition])
        if(direction == -1){
          amountBackward <- min(amountBackward, min(model[["I"]][(newPosition + 1):position]) - 1)
        }
        model_lp_initial <- model_lp_initial -
          (-
             #Choosing that time point
             log(sum(model[[target]]!=0)) -
             #choosing that that amount point of points to move
             log(
               amountBackward
               ) -
             #choosing the size of the move
             log(min(
               maxStep,
               (newPosition - 1)*(-direction == -1) +
                 (length(model[[target]]) - newPosition)*(-direction == 1)
             ) - minStep + 1) -
             #choosing direction
             log(1 + 1*(newPosition != length(model[[target]]) & newPosition != 1))
          )


        log_MH_ratio <- model_lp_proposed - model_lp_initial

        u <- runif(1, 0, 1)
        if(u < exp(log_MH_ratio)){
          model[["tracers"]][1:3,evalCol] <<- c(model[["tracers"]][1,evalCol] + 1/runs,
                                                model[["tracers"]][2,evalCol] + size^2/runs,
                                                model[["tracers"]][3,evalCol] +
                                                  mean((trueValue - model[[target]])^2)/runs)
          jump <- TRUE
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
        }
      },
      methods = list(
        reset = function () {}
      )
    )
  }}
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  control <- buildControl(epiModel, hyperParameters, "Beta", 1)
  output$addSampler(target = 'Beta', type = sampler_RW_TRACER, control = control)
  control <- buildControl(epiModel, hyperParameters, "Gamma", 2)
  output$addSampler(target = 'Gamma', type = sampler_RW_TRACER, control = control)
  control <- buildControl(epiModel, hyperParameters, "newI", 3)
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
  control <- buildControl(epiModel, hyperParameters, "Beta", 1)
  output$addSampler(target = 'Beta', type = sampler_RW_TRACER, control = control)
  control <- buildControl(epiModel, hyperParameters, "Gamma", 2)
  output$addSampler(target = 'Gamma', type = sampler_RW_TRACER, control = control)
  control <- buildControl(epiModel, hyperParameters, "newR", 3)
  output$addSampler(target = "newR",
                    type = stepSampler,
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
buildControl <- function(epiModel, hyperParameters, parameter, column){
  control <- hyperParameters[[parameter]]
  control$Column <- column
  if(identical(hyperParameters$`True Values`[[parameter]], NULL)){
    control$True <- rep(0, length(epiModel@Model[[parameter]]))
  }
  else{
    control$True <- hyperParameters$`True Values`[[parameter]]
  }
  return(control)
}
