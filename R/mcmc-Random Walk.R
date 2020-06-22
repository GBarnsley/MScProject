#' Generates a new candidate using the randomwalk method.
#' @export
randomWalk <- function(epiModel, priors){
  UseMethod("randomWalk")
}
#' Generates a new candidate using the randomwalk for SIR class.
#' Generates Beta, Gamma from normal centered around the past value of
#' Beta, Gamma and with standard deviation given in hyperparameters.
#' @export
randomWalk.SIR <- function(epiModel, hyperParameters){
  epiModel2 <- epiModel
  epiModel2@Beta <- rnorm(1,
                         mean = epiModel@Beta,
                         sd = hyperParameters$`Random Walk`$Beta)
  epiModel2@Gamma <- rnorm(1,
                          mean = epiModel@Gamma,
                          sd = hyperParameters$`Random Walk`$Gamma)
  if(epiModel2@Beta <= 0|epiModel2@Gamma <= 0){
    return(NA)
  }
  return(epiModel2)
}
#' Generates a new candidate using the randomwalk for iSIR class.
#' Generates Beta, Gamma from normal centered around the past value of
#' Beta, Gamma and with standard deviation given in hyperparameters.
#' @export
randomWalk.iSIR <- function(epiModel, hyperParameters){
  epiModel2 <- epiModel
  epiModel2@Beta <- rnorm(1,
                          mean = epiModel@Beta,
                          sd = hyperParameters$`Random Walk`$Beta)
  epiModel2@Gamma <- rnorm(1,
                           mean = epiModel@Gamma,
                           sd = hyperParameters$`Random Walk`$Gamma)
  return(simulate(epiModel2))
}
#' Generates a new candidate using the randomwalk for rSIR class.
#' Generates Beta, Gamma from normal centered around the past value of
#' Beta, Gamma and with standard deviation given in hyperparameters.
#' @export
randomWalk.rSIR <- function(epiModel, hyperParameters){
  return(randomWalk.iSIR(epiModel, hyperParameters))
}

#' Finds the log probability that the random walk generated the given candidate.
#' @export
logRandomWalk <- function(epiModel1, epiModel2){
  UseMethod("logRandomWalk")
}
#' Method for the iSIR class.
#' Doesn't find the density for Beta, Gamma since these will cancel
#' in the acceptance probability. Sums the log densities for newI if
#' Beta and Gamma came from the old model.
#' @export
logRandomWalk.iSIR <- function(epiModel1, epiModel2){
  ll <-
    0
  for(t in 1:min(sum(epiModel2@S!=0), length(epiModel2@newI))){
    ll <- ll +
      dbinom(epiModel2@newI[t],
             epiModel2@S[t],
             probGen(epiModel1@Beta*epiModel2@t.step*epiModel2@I[t]),
             log = TRUE)
  }
  return(
    ll
  )
}
#' Method for the iSIR class.
#' Doesn't find the density for Beta, Gamma since these will cancel
#' in the acceptance probability. Sums the log densities for newR if
#' Beta and Gamma came from the old model.
#' @export
logRandomWalk.rSIR <- function(epiModel1, epiModel2){
  ll <-
    0
  for(t in 1:min(sum(epiModel2@I!=0), length(epiModel2@newR))){
    ll <- ll +
      dbinom(epiModel2@newR[t],
             epiModel2@I[t],
             probGen(epiModel1@Gamma*epiModel2@t.step),
             log = TRUE)
  }
  return(
    ll
  )
}
#' Method for the SIR class.
#' Just returns 0 since its symmetric.
#' @export
logRandomWalk.SIR <- function(epiModel1, epiModel2){
  return(
    0
  )
}
