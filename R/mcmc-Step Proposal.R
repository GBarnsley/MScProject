#'
#' @export
stepProp <- function(epiModel, hyperParameters, i){
  UseMethod("stepProp")
}
#'
#' @export
stepProp.iSIR <- function(epiModel, hyperParameters, i){
  direction <- sample(c(-1,1), 1, replace = TRUE)
  if(direction == 1){
    lower <- sample(x = which(epiModel@newI!=0), size = 1, replace = TRUE)
    upper <- lower + sample(x = hyperParameters$`Step Proposal`$`Minimum Step Size`:hyperParameters$`Step Proposal`$`Maximum Step Size`,
                                              size = 1,
                                              replace = TRUE)
  }
  else{
    upper <- sample(x = which(epiModel@newI!=0), size = 1, replace = TRUE)
    lower <- upper - sample(x = hyperParameters$`Step Proposal`$`Minimum Step Size`:hyperParameters$`Step Proposal`$`Maximum Step Size`,
                            size = 1,
                            replace = TRUE)
  }
  if(lower <= 0|upper > length(epiModel@newI)){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  newEpiModel <- epiModel
  newEpiModel@newI[lower] <- epiModel@newI[lower] - direction
  newEpiModel@S[(lower+1):upper] <- epiModel@S[(lower+1):upper] + direction
  newEpiModel@I[(lower+1):upper] <- epiModel@I[(lower+1):upper] - direction
  if(any(newEpiModel@S[(lower+1):upper]<=0)|any(newEpiModel@I[(lower+1):upper]<=0)){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  newEpiModel@newI[upper] <- epiModel@newI[upper] + direction
  if(logPrior(newEpiModel, hyperParameters$Priors) +
     logLikelihood(newEpiModel) +
     -log(newEpiModel@newI[c(lower, upper)[direction/2 + 3/2]]) +
     -log(sum(newEpiModel@newI != 0)) -
     logPrior(epiModel, hyperParameters$Priors) -
     logLikelihood(epiModel) -
     -log(epiModel@newI[c(upper, lower)[direction/2 + 3/2]]) -
     -log(sum(epiModel@newI != 0)) >=
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
#' @export
stepProp.rSIR <- function(epiModel, hyperParameters, i){
  direction <- sample(c(-1,1), 1, replace = TRUE)
  if(direction == 1){
    lower <- sample(x = which(epiModel@newR!=0), size = 1, replace = TRUE)
    upper <- lower + sample(x = hyperParameters$`Step Proposal`$`Minimum Step Size`:hyperParameters$`Step Proposal`$`Maximum Step Size`,
                            size = 1,
                            replace = TRUE)
  }
  else{
    upper <- sample(x = which(epiModel@newR!=0), size = 1, replace = TRUE)
    lower <- upper - sample(x = hyperParameters$`Step Proposal`$`Minimum Step Size`:hyperParameters$`Step Proposal`$`Maximum Step Size`,
                            size = 1,
                            replace = TRUE)
  }
  if(lower <= 0|upper > length(epiModel@newR)){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  newEpiModel <- epiModel
  newEpiModel@newR[lower] <- epiModel@newR[lower] - direction
  newEpiModel@I[(lower+1):upper] <- epiModel@I[(lower+1):upper] + direction
  newEpiModel@R[(lower+1):upper] <- epiModel@R[(lower+1):upper] - direction
  if(any(newEpiModel@I[(lower+1):upper]<=0)|any(newEpiModel@R[(lower+1):upper]<=0)){
    epiModel@attempts <- epiModel@attempts + 1
    return(epiModel)
  }
  newEpiModel@newR[upper] <- epiModel@newR[upper] + direction
  if(logPrior(newEpiModel, hyperParameters$Priors) +
     logLikelihood(newEpiModel) +
     -log(newEpiModel@newR[c(lower, upper)[direction/2 + 3/2]]) +
     -log(sum(newEpiModel@newR != 0)) -
     logPrior(epiModel, hyperParameters$Priors) -
       logLikelihood(epiModel) -
       -log(epiModel@newR[c(upper, lower)[direction/2 + 3/2]]) -
       -log(sum(epiModel@newR != 0))>=
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
