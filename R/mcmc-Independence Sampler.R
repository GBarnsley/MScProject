#' Generates a candidate sample from a predetermined fixed distribution
#' @param epiModel epidemic model
#' @param hyperParameters list of hyperparameters, used in distribution
#' @export
independenceSampler <- function(epiModel, hyperParameters){
  UseMethod("independenceSampler")
}
#' Independence sampler for the SIR class.
#' First generates beta and gamma from a normal distribution centred around the mean
#' value of their priors. Then uses those values to generate a newI or newR using
#' the standard Binomial distribution from the SIR model.
#' @export
independenceSampler.SIR <- function(epiModel, hyperParameters){
  epiModel@Beta <- rnorm(1,
                         mean = hyperParameters$`Independence Sampler`$Beta$Mean,
                         sd = hyperParameters$`Independence Sampler`$Beta$SD)
  epiModel@Gamma <- rnorm(1,
                          mean = hyperParameters$`Independence Sampler`$Gamma$Mean,
                          sd = hyperParameters$`Independence Sampler`$Gamma$SD)
  if(epiModel@Beta <= 0|epiModel@Gamma <= 0){
    return(NA)
  }
  return(epiModel)
}
#' Independence sampler for the iSIR class.
#' First generates beta and gamma from a normal distribution centred around the mean
#' value of their priors. Then uses those values to generate a newI or newR using
#' the standard Binomial distribution from the SIR model.
#' @export
independenceSampler.iSIR <- function(epiModel, hyperParameters){
  epiModel@Beta <- rnorm(1,
                         mean = hyperParameters$`Independence Sampler`$Beta$Mean,
                         sd = hyperParameters$`Independence Sampler`$Beta$SD)
  epiModel@Gamma <- rnorm(1,
                          mean = hyperParameters$`Independence Sampler`$Gamma$Mean,
                          sd = hyperParameters$`Independence Sampler`$Gamma$SD)
  return(simulate(epiModel))
}
#' Independence sampler for the rSIR class.
#' First generates beta and gamma from a normal distribution centred around the mean
#' value of their priors. Then uses those values to generate a newI or newR using
#' the standard Binomial distribution from the SIR model.
#' @export
independenceSampler.rSIR <- function(epiModel, hyperParameters){
  return(simulate,iSIR(epiModel, hyperParameters))
}
#' Finds the log probability of that the indepdence sampler genreated
#' that sample.
#' @export
logIndependenceSampler <- function(epiModel, hyperParameters){
  UseMethod("logIndependenceSampler")
}
#' Method for the iSIR class.
#' finds the log density for the values of Beta and Gamma being
#' generated and then that the binomial chain generated the values of
#' newI.
#' @export
logIndependenceSampler.iSIR <- function(epiModel, hyperParameters){
  ll <-
    dnorm(epiModel@Beta,
          mean = hyperParameters$Beta$Mean,
          sd = hyperParameters$Beta$SD,
          log = TRUE) +
    dnorm(epiModel@Gamma,
          mean = hyperParameters$Gamma$Mean,
          sd = hyperParameters$Gamma$SD,
          log = TRUE)
  for(t in 1:min(sum(epiModel@S != 0), length(epiModel@newI))){
    ll <- ll +
      dbinom(epiModel@newI[t],
             epiModel@S[t],
             probGen(epiModel@Beta*epiModel@t.step*epiModel@I[t]),
             log = TRUE)
  }
  return(
    ll
  )
}
#' Method for the iSIR class.
#' finds the log density for the values of Beta and Gamma being
#' generated and then that the binomial chain generated the values of
#' newR.
#' @export
logIndependenceSampler.rSIR <- function(epiModel, hyperParameters){
  ll <-
    dnorm(epiModel@Beta,
          mean = hyperParameters$Beta$Mean,
          sd = hyperParameters$Beta$SD,
          log = TRUE) +
    dnorm(epiModel@Gamma,
          mean = hyperParameters$Gamma$Mean,
          sd = hyperParameters$Gamma$SD,
          log = TRUE)
  for(t in 1:min(sum(epiModel@I!=0), length(epiModel@newR))){
    ll <- ll +
       dbinom(epiModel@newR[t],
              epiModel@I[t],
              probGen(epiModel@Gamma*epiModel@t.step),
              log = TRUE)
  }
  return(
    ll
  )
}
#' Method for the SIR class.
#' finds the log density for the values of Beta and Gamma being
#' generated.
#' @export
logIndependenceSampler.SIR <- function(epiModel, hyperParameters){
  return(
    dnorm(epiModel@Beta,
          mean = hyperParameters$Beta$Mean,
          sd = hyperParameters$Beta$SD,
          log = TRUE) +
    dnorm(epiModel@Gamma,
          mean = hyperParameters$Gamma$Mean,
          sd = hyperParameters$Gamma$SD,
          log = TRUE)
  )
}
