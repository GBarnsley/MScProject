#'
#' @export
<<<<<<< HEAD
stepProp <- function(epiModel, hyperParameters, i){
=======
stepProp <- function(epiModel, hyperParameters){
>>>>>>> 1c714e3e41b618fd584d43c4f5ae32a0417f471a
  UseMethod("stepProp")
}
#'
#' @export
<<<<<<< HEAD
stepProp.iSIR <- function(epiModel, hyperParameters, i){
=======
stepProp.iSIR <- function(epiModel, hyperParameters){
>>>>>>> 1c714e3e41b618fd584d43c4f5ae32a0417f471a
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
<<<<<<< HEAD
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
=======
    return(NA)
  }
  epiModel@newI[lower] <- epiModel@newI[lower] - direction
  epiModel@S[(lower+1):upper] <- epiModel@S[(lower+1):upper] + direction
  epiModel@I[(lower+1):upper] <- epiModel@I[(lower+1):upper] - direction
  if(any(epiModel@S[(lower+1):upper]<=0)|any(epiModel@I[(lower+1):upper]<=0)){
    return(NA)
  }
  epiModel@newI[upper] <- epiModel@newI[upper] + direction
  return(epiModel)
}
#'
#' @export
stepProp.rSIR <- function(epiModel, hyperParameters){
>>>>>>> 1c714e3e41b618fd584d43c4f5ae32a0417f471a
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
<<<<<<< HEAD
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
=======
    return(NA)
  }
  epiModel@newR[lower] <- epiModel@newR[lower] - direction
  epiModel@I[(lower+1):upper] <- epiModel@I[(lower+1):upper] + direction
  epiModel@R[(lower+1):upper] <- epiModel@R[(lower+1):upper] - direction
  if(any(epiModel@I[(lower+1):upper]<=0)|any(epiModel@R[(lower+1):upper]<=0)){
    return(NA)
  }
  epiModel@newR[upper] <- epiModel@newR[upper] + direction
  return(epiModel)
}
#'
#' @export
logStep <- function(epiModel1, epiModel2, hyperParameters){
  UseMethod("logStep")
}
#'
#' @export
logStep.iSIR <- function(epiModel1, epiModel2){
  position <- which(epiModel2@newI - epiModel1@newI == -1)
  return(
    -log(epiModel1@newI[position])
  )
}
#'
#' @export
logStep.rSIR <- function(epiModel1, epiModel2){
  position <- which(epiModel2@newR - epiModel1@newR == -1)
  return(
    -log(epiModel1@newR[position])
  )
>>>>>>> 1c714e3e41b618fd584d43c4f5ae32a0417f471a
}
