#' Generates a candidate sample from a predetermined fixed distribution
#' @param epiModel epidemic model
#' @param hyperParameters list of hyperparameters, used in distribution
#' @export
independenceSampler <- function(epiModel, hyperParameters, i){
  UseMethod("independenceSampler")
}
#' Independence sampler for the SIR class.
#' First generates beta and gamma from a normal distribution centred around the mean
#' value of their priors.
#' @export
independenceSampler.SIR <- function(epiModel, hyperParameters, i){
  return(
    independenceSamplerGamma(
      independenceSamplerBeta(
        epiModel, hyperParameters
      ),
      hyperParameters
    )
  )
}
#'
#'@export
independenceSamplerBeta <- function(epiModel, hyperParameters){
  newEpiModel <- epiModel
  newEpiModel@Beta <- rnorm(1,
                            mean = hyperParameters$`Independence Sampler`$Beta$Mean,
                            sd = hyperParameters$`Independence Sampler`$Beta$SD)
  if(newEpiModel@Beta <= 0){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  else if(logPrior(newEpiModel, hyperParameters$Priors) +
          logLikelihood(newEpiModel) +
          dnorm(epiModel@Beta,
                mean = hyperParameters$`Independence Sampler`$Beta$Mean,
                sd = hyperParameters$`Independence Sampler`$Beta$SD,
                log = TRUE) -
          logPrior(epiModel, hyperParameters$Priors) -
          logLikelihood(epiModel) -
          dnorm(newEpiModel@Beta,
                mean = hyperParameters$`Independence Sampler`$Beta$Mean,
                sd = hyperParameters$`Independence Sampler`$Beta$SD,
                log = TRUE) >=
          log(runif(1))){
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
independenceSamplerGamma <- function(epiModel, hyperParameters){
  newEpiModel <- epiModel
  newEpiModel@Gamma <- rnorm(1,
                            mean = hyperParameters$`Independence Sampler`$Gamma$Mean,
                            sd = hyperParameters$`Independence Sampler`$Gamma$SD)
  if(newEpiModel@Gamma <= 0){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  else if(logPrior(newEpiModel, hyperParameters$Priors) +
          logLikelihood(newEpiModel) +
          dnorm(epiModel@Gamma,
                mean = hyperParameters$`Independence Sampler`$Gamma$Mean,
                sd = hyperParameters$`Independence Sampler`$Gamma$SD,
                log = TRUE) -
          logPrior(epiModel, hyperParameters$Priors) -
          logLikelihood(epiModel) -
          dnorm(newEpiModel@Beta,
                mean = hyperParameters$`Independence Sampler`$Gamma$Mean,
                sd = hyperParameters$`Independence Sampler`$Gamma$SD,
                log = TRUE) >=
          log(runif(1))){
    newEpiModel@acceptance <- epiModel@acceptance + 1
    return(newEpiModel)
  }
  else{
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
}
