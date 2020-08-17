#' Plots some useful graphs of the SIR Model
#'
#' Plots the populations over time
#' Requires ggplot
#' Won't work currently, needs rewriting
#'
#' @param epiModel An epidemic model of SIR class
#' @return ggplot object
#' @export
plot.SIR <- function(epiModel){
  n <- length(epiModel@S)
  temp <- data.frame(Population = rep(0,3*n), State = rep(NA,3*n), Time = rep(0,3*n))
  for(i in 1:n){
    temp[3*(i-1) + 1,] <- c(epiModel@S[i], "Susceptible", i)
    temp[3*(i-1) + 2,] <- c(epiModel@I[i], "Infectious", i)
    temp[3*(i-1) + 3,] <- c(epiModel@R[i], "Recovered", i)
  }
  temp$Population <- as.numeric(temp$Population)
  #temp$State <- as.factor(temp$State)
  temp$Time <- as.numeric(temp$Time)
  print(
    ggplot2::ggplot(data = temp,
                    ggplot2::aes(x = Time,
                        y = Population,
                        colour = State))+
      ggplot2::geom_line(size = 1, alpha = 0.8) +
      ggplot2::scale_colour_manual(values = c("Susceptible" = "green",
                                              "Infectious" = "red",
                                              "Recovered" = "black"))
  )
}
#' sa;dlsal;df
#' @export
plotSamples <- function(epiModel, timePoints){
  UseMethod("plotSamples")
}
#' Plots some useful graphs of the bayesian-SIR Model
#'
#' Plots the populations over time
#' Requires ggplot
#'
#' @param samples An object containing many epidemic model of SIR class
#' @return ggplot object
#' @export
plotSamples.SIR <- function(epiModel, timePoints){
  print(ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(
            x = 1:nrow(epiModel@Samples),
            y = epiModel@Samples[,'Beta']
          )) +
            ggplot2::labs(x = "Index", y = "Beta"))
  print(ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(
            x = 1:nrow(epiModel@Samples),
            y = epiModel@Samples[,'Gamma']
          )) +
            ggplot2::labs(x = "Index", y = "Gamma"))
}
#' Plots some useful graphs of the bayesian-iSIR Model
#'
#' Plots the populations over time
#' Requires ggplot
#'
#' @param samples An object containing many epidemic model of SIR class
#' @return ggplot object
#' @export
plotSamples.iSIR <- function(epiModel, timePoints = NA){
  print(ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(
            x = 1:nrow(epiModel@Samples),
            y = epiModel@Samples[,'Beta']
          )) +
            ggplot2::labs(x = "Index", y = "Beta"))
  print(ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(
            x = 1:nrow(epiModel@Samples),
            y = epiModel@Samples[,'Gamma']
          )) +
            ggplot2::labs(x = "Index", y = "Gamma"))
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(epiModel)
  }
  for(time in timePoints){
    print(ggplot2::ggplot() +
            ggplot2::geom_line(ggplot2::aes(
              x = 1:nrow(epiModel@Samples),
              y = epiModel@Samples[,paste0("newI[",time,"]")]
            )) +
              ggplot2::labs(x = "Index", y = paste0("new I at time ",time)))
  }
}
#' Plots some useful graphs of the bayesian-rSIR Model
#'
#' Plots the populations over time
#' Requires ggplot
#'
#' @param samples An object containing many epidemic model of SIR class
#' @return ggplot object
#' @export
plotSamples.rSIR <- function(epiModel, timePoints = NA){
  print(ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(
            x = 1:nrow(epiModel@Samples),
            y = epiModel@Samples[,'Beta']
          )) +
            ggplot2::labs(x = "Index", y = "Beta"))
  print(ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(
            x = 1:nrow(epiModel@Samples),
            y = epiModel@Samples[,'Gamma']
          )) +
            ggplot2::labs(x = "Index", y = "Gamma"))
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(epiModel)
  }
  for(time in timePoints){
    print(ggplot2::ggplot() +
            ggplot2::geom_line(ggplot2::aes(
              x = 1:nrow(epiModel@Samples),
              y = epiModel@Samples[,paste0("newR[",time,"]")]
            )) +
              ggplot2::labs(x = "Index", y = paste0("new R at time ",time)))
  }
}
