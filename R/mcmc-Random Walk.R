#' Generates a new candidate using the randomwalk method.
#' @export
randomWalk <- function(epiModel, hyperParameters, i){
  UseMethod("randomWalk")
}
#' Generates a new candidate using the randomwalk for SIR class.
#' Generates Beta, Gamma from normal centered around the past value of
#' Beta, Gamma and with standard deviation given in hyperparameters.
#' @export
randomWalk.SIR <- function(epiModel, hyperParameters, i){
  if(i%%2 == 0){
    randomWalkBeta(
      epiModel,
      hyperParameters
    )
  }
  else{
    randomWalkGamma(
      epiModel,
      hyperParameters
    )
  }
}
#'
#'@export
randomWalkBeta <- function(epiModel, hyperParameters){
  newEpiModel <- epiModel
  newEpiModel@Beta <- rnorm(1,
                            mean = epiModel@Beta,
                            sd = hyperParameters$`Random Walk`$Beta)
  if(newEpiModel@Beta <= 0){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  else if(logPrior(newEpiModel, hyperParameters$Priors) +
          logLikelihood(newEpiModel) -
          logPrior(epiModel, hyperParameters$Priors) -
          logLikelihood(epiModel) >=
          log(runif(1))){
    newEpiModel@attempts <- epiModel@attempts + 1
    newEpiModel@acceptance <- epiModel@acceptance + 1
    return(newEpiModel)
  }
  else{
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
}
#'
#'@export
randomWalkGamma <- function(epiModel, hyperParameters){
  newEpiModel <- epiModel
  newEpiModel@Gamma <- rnorm(1,
                            mean = epiModel@Gamma,
                            sd = hyperParameters$`Random Walk`$Gamma)
  if(newEpiModel@Gamma <= 0){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  else if(logPrior(newEpiModel, hyperParameters$Priors) +
          logLikelihood(newEpiModel) -
          logPrior(epiModel, hyperParameters$Priors) -
          logLikelihood(epiModel) >=
          log(runif(1))){
    newEpiModel@attempts <- epiModel@attempts + 1
    newEpiModel@acceptance <- epiModel@acceptance + 1
    return(newEpiModel)
  }
  else{
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
}
