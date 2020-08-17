#'
#'@export
stepSampler_setup <- function(model, mvSaved, target, control) {
  maxStep <- as.integer(control$TMax)
  maxChange <- as.integer(control$DeltaMax)
  calcNodes <- model$getDependencies(target)
  runs <- as.integer(control$R)
}
#'
#'@export
stepSampler_run <- function() {
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
    amount <- rcat(n = 1, prob = rep(1, amounts))

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
    if(u < exp(log_MH_ratio)){
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
  }
}
#'
#'@export
stepSampler_run_bounded <- function() {
  model_lp_initial <- getLogProb(model, calcNodes)
  positions <- which(model[[target]]!=0)
  for(i in 1:runs){
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

    sizes <- min(maxStep, (position - 1)*(direction == -1) +
                   (length(model[[target]]) - position)*(direction == 1))
    size <- rcat(n = 1, prob = rep(1, sizes))

    newPosition <- position + direction*size

    amounts <- min(maxChange, model[[target]][position])
    if(direction == 1){
      amounts <- min(amounts, min(model[["I"]][(position + 1):newPosition]) - 1)
    }
    if(amounts > 0){
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
      amountBackward <- min(maxChange, model[[target]][newPosition])
      if(direction == -1){
        amountBackward <- min(amountBackward, min(model[["I"]][(newPosition + 1):position]) - 1)
      }
      if(amountBackward > 0){
        sampler_lp_initial <-
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
             )) -
             #choosing direction
             log(1 + 1*(newPosition != length(model[[target]]) & newPosition != 1))
          )
        log_MH_ratio <- (model_lp_proposed - sampler_lp_proposed) - (model_lp_initial - sampler_lp_initial)

        u <- runif(1, 0, 1)
        if(u < exp(log_MH_ratio)){
          model_lp_initial <- model_lp_proposed
          positions <- which(model[[target]]!=0)
          jump <- TRUE
        }
        else{
          jump <- FALSE
        }
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
  }
}
