#'
#'@export
buildMCMCInternal <- function(epiModel, hyperParameters){
  stepSampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
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
      positions <- which(model[[target]]!=0)
      for(i in 1:runs){
        #Choose the original position
        pos <- rcat(n = 1, prob = rep(1/length(positions), length(positions)))
        position <- positions[pos]
        #Choosing the direction of the move
        if(position == length(model[[target]])){
          #if the position is the last time point we can only go backward
          directions <- c(-1)
          direction <- -1
        }
        else if(position == 1){
          #if it is the first time point we have to go forward
          directions <- c(1)
          direction <- 1
        }
        else{
          #in all other cases we choose with uniform probability either direction
          directions <- c(-1,1)
          pos <- rcat(n = 1, prob = c(0.5,0.5))
          direction <- directions[pos]
        }
        #Choosing the number of places to move
        sizes <- min(maxStep, (position - 1)*(direction == -1) +
                       (length(model[[target]]) - position)*(direction == 1))
        #prevents moving further than the time frame of the epidemic
        size <- rcat(n = 1, prob = rep(1/sizes, sizes))

        #the index of the position we are moving to
        newPosition <- position + direction*size

        #choosing the number of points to move
        amounts <- min(maxChange, model[[target]][position])
        #stops us from moving more than the number of points at that time
        amount <- rcat(n = 1, prob = rep(1/amounts, amounts))

        sampler_lp_proposed <-
          (-
             #Choosing that time point
             log(length(positions)) -
             #choosing that that amount point of points to move
             log(amounts) -
             #Choosing the size of the move
             log(sizes) -
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
             log(
               min(maxChange, model[[target]][newPosition])
             ) -
             #choosing the size of the move
             log(min(
               maxStep,
               (newPosition - 1)*(-direction == -1) +
                 (length(model[[target]]) - newPosition)*(-direction == 1)
             )) -
             #choosing direction
             log(1 + 1*(newPosition != length(model[[target]]) & newPosition != 1))
          )
        log_MH_ratio <- (model_lp_proposed - sampler_lp_proposed) - (model_lp_initial - sampler_lp_initial)

        u <- runif(1, 0, 1)
        model[["tracers"]][1,evalCol] <<- model[["tracers"]][1,evalCol] + min(exp(log_MH_ratio),1)/runs
        if(u < exp(log_MH_ratio)){
          model[["tracers"]][2,evalCol] <<- model[["tracers"]][2,evalCol] + (size*amount)^2/runs
          model_lp_initial <- model_lp_proposed
          positions <- which(model[[target]]!=0)
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
        model[["tracers"]][3,evalCol] <<- model[["tracers"]][3,evalCol] +
          mean((trueValue - model[[target]])^2)/runs
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
        maxStep <- as.integer(control$`Maximum Step Size`)
        maxChange <- as.integer(control$`Maximum Change`)
        runs <- as.integer(control$`Runs per Random Walk`)
        calcNodes <- model$getDependencies(target)
        evalCol <- as.integer(control$Column)
        trueValue <- control$True
      },
      ## the run function
      run = function() {
        model_lp_initial <- getLogProb(model, calcNodes)
        pointsToMove <- model[[target]]!=0
        canMoveBackward <- rep(TRUE, length(model[[target]]))
        canMoveBackward[1] <- FALSE
        canMoveForward <- rep(FALSE, length(model[[target]]))
        for(i in 1:(length(model[[target]])-1)){
          canMoveForward[i] <- model[["I"]][i+1] > 1
        }
        positions <- which(pointsToMove * (canMoveBackward + canMoveForward - canMoveBackward*canMoveForward))
        jump <- FALSE

        for(i in 1:runs){

          #Choose the original position
          pos <- rcat(n = 1, prob = rep(1, length(positions)))
          position <- positions[pos]

          #Choosing the direction of the move
          directions <- canMoveBackward[position] + canMoveForward[position]
          direction <- c(-1,1)[rcat(n=1, prob = c(canMoveBackward[position], canMoveForward[position]))]


          #Choosing the number of places to move
          if(direction == 1){
            sizes <- min(maxStep, length(model[[target]]) - position)
            #now we check to make sure we can't propose a step where the amount <= 0
            minAmounts <- rep(0, sizes)
            for(i in 1:sizes){
              minAmounts[i] <- min(model[["I"]][(position + 1):(position + i)]) - 1
            }
            #if there is a value where amount <= 0
            if(min(minAmounts) == 0){
              sizes <- min(which(minAmounts == 0)) - 1
              #then the step size before this occurs is our new maximum
            }
          }
          else{
            sizes <- min(maxStep, position - 1)
          }
          size <- rcat(n = 1, prob = rep(1, sizes))

          #the index of the position we are moving to
          newPosition <- position + direction*size

          #choosing the number of points to move
          amounts <- min(maxChange, model[[target]][position])
          #stops us from moving more than the number of points at that time
          if(direction == 1){
            amounts <- min(amounts, minAmounts[size])
            #also bounded by the minimum amounts for that size, which has already been calculated
          }
          amount <- rcat(n = 1, prob = rep(1, amounts))

          #the probability of generating this proposal
          sampler_lp_proposed <-
            (-
               #Choosing that time point
               log(length(positions)) -
               #choosing that that amount point of points to move
               log(amounts) -
               #Choosing the size of the move
               log(sizes) -
               #choosing direction
               log(directions)
            )

          #generating some values for the reverse probability
          #some of these are easier to find when we haven't updated the values
          positionsRev <- length(positions) + (model[[target]][newPosition] == 0) - (amount == model[[target]][position])

          model[[target]][position] <<- model[[target]][position] - amount
          model[[target]][newPosition] <<- model[[target]][newPosition] + amount
          model_lp_proposed <- calculate(model, calcNodes)

          #Calculating the rest of the values for the reverse probability
          if(newPosition == length(model[[target]])){
            directionsRev <- 1
            #if we are at the last timepoint we can only go backward
          }
          else{
            directionsRev <-
              (newPosition != 1) +
              #backwards move always possible unless at the start
              (model[["I"]][newPosition + 1] > 1)
              #is a forward move possible
          }
          amountsRev <- min(maxChange, model[[target]][newPosition])
          if(direction == -1){
            amountsRev <- min(amountsRev, min(model[["I"]][(newPosition+1):(newPosition + size)]) - 1)
          }
          if(direction == -1){
            sizesRev <- min(maxStep, length(model[[target]]) - newPosition)
            if(sizesRev > size){
              #we know that it must be possible for the size of the proposed jump
              minAmounts <- rep(Inf, sizesRev)
              for(i in (size+1):sizesRev){
                minAmounts[i] <- min(model[["I"]][(newPosition + size + 1):(newPosition + i)]) - 1
              }
              #if there is a value where amount <= 0
              if(min(minAmounts) == 0){
                sizesRev <- min(which(minAmounts == 0)) - 1
                #then the step size before this occurs is our new maximum
              }
            }
          }
          else{
            sizesRev <- min(maxStep, newPosition - 1)
          }

          #probability of generating the initial values from the proposed state
          sampler_lp_initial <-
            (-
               #Choosing that time point
               log(positionsRev) -
               #choosing that that amount point of points to move
               log(amountsRev) -
               #choosing the size of the move
               log(sizesRev) -
               #choosing direction
              log(directionsRev)
          )
          log_MH_ratio <- (model_lp_proposed - sampler_lp_proposed) - (model_lp_initial - sampler_lp_initial)

          model[["tracers"]][1,evalCol] <<- model[["tracers"]][1,evalCol] + min(exp(log_MH_ratio),1)/runs
          u <- runif(1, 0, 1)
          if(u < exp(log_MH_ratio)){
            model[["tracers"]][2,evalCol] <<- model[["tracers"]][2,evalCol] + (size*amount)^2/runs
            model_lp_initial <- model_lp_proposed
            pointsToMove <- model[[target]]!=0
            for(i in 1:(length(model[[target]])-1)){
              canMoveForward[i] <- model[["I"]][i+1] > 1
            }
            positions <- which(pointsToMove * (canMoveBackward + canMoveForward - canMoveBackward*canMoveForward))

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
          model[["tracers"]][3,evalCol] <<- model[["tracers"]][3,evalCol] +
            mean((trueValue - model[[target]])^2)/runs
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
