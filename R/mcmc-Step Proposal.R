#'
#' @export
stepProp <- function(epiModel, hyperParameters){
  UseMethod("stepProp")
}
#'
#' @export
stepProp.iSIR <- function(epiModel, hyperParameters){
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
}
