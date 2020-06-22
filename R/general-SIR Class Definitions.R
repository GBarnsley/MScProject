SIRclass <- setClass(
  "SIR",
  slots = c(
    newI = "ANY",
    newR = "ANY",
    N = "numeric",
    Beta = "ANY",
    Gamma = "ANY",
    S = "ANY",
    I = "ANY",
    R = "ANY",
    t.step = "numeric"
    )
)
newISIRclass <- setClass(
  "iSIR",
  slots = c(
    newI = "ANY",
    newR = "ANY",
    N = "numeric",
    Beta = "ANY",
    Gamma = "ANY",
    S = "ANY",
    I = "ANY",
    R = "ANY",
    t.step = "numeric"
  ),
  contains = "SIR"
)
newRSIRclass <- setClass(
  "rSIR",
  slots = c(
    newI = "ANY",
    newR = "ANY",
    N = "numeric",
    Beta = "ANY",
    Gamma = "ANY",
    S = "ANY",
    I = "ANY",
    R = "ANY",
    t.step = "numeric"
  ),
  contains = "SIR"
)
#' Creates an object of class SIR, iSIR or rSIR
#'
#' Takes any combination of the values in the SIR model and uses dataGen to generate
#' unspecified data where it can. If S, I and newI but R and newR are then it will
#' generate an object of class iSIR that inherits from the SIR class. This class functions
#' the same for most functions other than in the metropolis hasting algorithm where it treats
#' newI as a parameter to generate. if R, I and newR are missing then it will make an object of
#' class rSIR which does the same but treats newR as a parameter.
#'
#' @param S Time series of the Susceptible Population
#' @param I Time series of the Infectious Population
#' @param R Time series of the Recovered Population
#' @param newI Time series of new infections
#' @param newR Time series of new recoveries
#' @param N The size of the population
#' @param Beta The infectious weight of one person in I
#' @param Gamma The recovery rate
#' @param t.step The amount of time in each step of the model, actual value depends on the units of Beta and Gamma
#' @return Object of class SIR
#' @export
SIR <- function(S = NULL,
                I = NULL,
                R = NULL,
                newI = NULL,
                newR = NULL,
                N = NULL,
                Beta = NULL,
                Gamma = NULL,
                t.step = 1){
  if(is.null(S)&is.null(I)&is.null(newI)&(!is.null(R)|!is.null(newR))){
    return(dataGen(newISIRclass(S = S,
                            I = I,
                            R = R,
                            newI = newI,
                            newR = newR,
                            N = N,
                            Beta = Beta,
                            Gamma = Gamma,
                            t.step = t.step)))
  }
  else if(is.null(R)&is.null(I)&is.null(newR)&(!is.null(S)|!is.null(newI))){
    return(dataGen(newRSIRclass(S = S,
                                I = I,
                                R = R,
                                newI = newI,
                                newR = newR,
                                N = N,
                                Beta = Beta,
                                Gamma = Gamma,
                                t.step = t.step)))
  }
  else{
    return(dataGen(SIRclass(S = S,
                            I = I,
                            R = R,
                            newI = newI,
                            newR = newR,
                            N = N,
                            Beta = Beta,
                            Gamma = Gamma,
                            t.step = t.step)))

  }
}
