setClass("samples")
#' Plots some useful graphs of the SIR Model
#'
#' Plots the populations over time
#' Requires ggplot
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
plotSamples <- function(x, timePoints){
  UseMethod("plotSamples", x)
}
#' Plots some useful graphs of the bayesian-SIR Model
#'
#' Plots the populations over time
#' Requires ggplot
#'
#' @param samples An object containing many epidemic model of SIR class
#' @return ggplot object
#' @export
plotSamples.SIR <- function(x, timePoints){
  data <- data.frame(Index = NA, Beta = NA, Gamma = NA)
  for(i in 1:length(x)){
    data[i,] <- c(i, x[[i]]@Beta, x[[i]]@Gamma)
  }
  print(ggplot2::ggplot(data = data,
                        ggplot2::aes(
                          x = Index,
                          y = Beta
                        )) +
          ggplot2::geom_line()
        )
  print(ggplot2::ggplot(data = data,
                              ggplot2::aes(
                                x = Index,
                                y = Gamma
                              )) +
          ggplot2::geom_line()
        )
}
#' Plots some useful graphs of the bayesian-iSIR Model
#'
#' Plots the populations over time
#' Requires ggplot
#'
#' @param samples An object containing many epidemic model of SIR class
#' @return ggplot object
#' @export
plotSamples.iSIR <- function(x, timePoints = NA){
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(x)
  }
  data <- matrix(NA, nrow = length(x), ncol = 1 + 2 + length(timePoints))
  for(i in 1:nrow(data)){
    data[i,1] <- i
    data[i,2] <- x[[i]]@Beta
    data[i,3] <- x[[i]]@Gamma
    for(j in 1:length(timePoints)){
      data[i,3+j] <- x[[i]]@newI[timePoints[j]]
    }
  }
  print(ggplot2::ggplot() +
          ggplot2::geom_line(
            ggplot2::aes(
              x = data[,1],
              y = data[,2]
            )
          ) +
          ggplot2::labs(
            x = "Index",
            y = "Beta"
          )
  )
  print(ggplot2::ggplot() +
          ggplot2::geom_line(
            ggplot2::aes(
              x = data[,1],
              y = data[,3]
            )
          ) +
          ggplot2::labs(
            x = "Index",
            y = "Gamma"
          )
  )
  for(i in 1:length(timePoints)){
    print(ggplot2::ggplot() +
            ggplot2::geom_line(
              ggplot2::aes(
                x = data[,1],
                y = data[,3+i]
              )
            ) +
            ggplot2::labs(
              x = "Index",
              y = paste("New I at time",
                        timePoints[i])
            )
    )
  }
  #data <- data.frame(Index = c(1:length(x[[i]]@R), rep(NA, length(x)*length(x[[1]]@S)*2)),
  #                   Population = c(x[[i]]@R, rep(NA, length(x)*length(x[[1]]@S)*2)),
  #                   State = c(rep("R", length(x[[i]]@R)), rep(NA, length(x)*length(x[[1]]@S)*2)),
  #                   ModelIndex = c(rep(0, length(x[[i]]@R)), rep(NA, length(x)*length(x[[1]]@S)*2)),
  #                   stringsAsFactors = FALSE)
  #for(i in 1:length(x)){
  #  for(j in 1:length(x[[i]]@S)){
  #    for(k in 1:2){
  #      if(k == 1){
  #        data[(i-1)*length(x[[i]]@S)*2 + j + k*length(x[[i]]@S),] <- c(j, x[[i]]@S[j], "S", i)
  #      }
  #      else{
  #        data[(i-1)*length(x[[i]]@S)*2 + j + k*length(x[[i]]@S),] <- c(j, x[[i]]@I[j], "I", i)
  #      }
  #    }
  #  }
  #}
  #data$Population <- as.numeric(data$Population)
  #data$Index <- as.numeric(data$Index)
  #data$State <- as.factor(data$State)
  #data$ModelIndex <- as.factor(data$ModelIndex)
  #print(
  #  ggplot2::ggplot(data) +
  #    ggplot2::geom_line(
  #      ggplot2::aes(x = Index,
  #                   y = Population,
  #                   colour = State,
  #                   group = ModelIndex:State),
  #      alpha = 0.2) +
  #    ggplot2::scale_colour_manual(values = c("S" = "green",
  #                                            "I" = "red",
  #                                            "R" = "black"))
  #)
  #plot(summarise(x, average))
}
#' Plots some useful graphs of the bayesian-rSIR Model
#'
#' Plots the populations over time
#' Requires ggplot
#'
#' @param samples An object containing many epidemic model of SIR class
#' @return ggplot object
#' @export
plotSamples.rSIR <- function(x, timePoints = NA){
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(x)
  }
  data <- matrix(NA, nrow = length(x), ncol = 1 + 2 + length(timePoints))
  for(i in 1:nrow(data)){
    data[i,1] <- i
    data[i,2] <- x[[i]]@Beta
    data[i,3] <- x[[i]]@Gamma
    for(j in 1:length(timePoints)){
      data[i,3+j] <- x[[i]]@newR[timePoints[j]]
    }
  }
  print(ggplot2::ggplot() +
          ggplot2::geom_line(
            ggplot2::aes(
              x = data[,1],
              y = data[,2]
            )
          ) +
          ggplot2::labs(
            x = "Index",
            y = "Beta"
          )
  )
  print(ggplot2::ggplot() +
          ggplot2::geom_line(
            ggplot2::aes(
              x = data[,1],
              y = data[,3]
            )
          ) +
          ggplot2::labs(
            x = "Index",
            y = "Gamma"
          )
  )
  for(i in 1:length(timePoints)){
    print(ggplot2::ggplot() +
            ggplot2::geom_line(
              ggplot2::aes(
                x = data[,1],
                y = data[,3+i]
              )
            ) +
            ggplot2::labs(
              x = "Index",
              y = paste("New R at time",
                        timePoints[i])
            )
    )
  }
  #data <- data.frame(Index = 1:length(x[[i]]@S),
  #                   Population = x[[i]]@S,
  #                   State = rep("S", length(x[[i]]@S)),
  #                   ModelIndex = rep(0, length(x[[i]]@S)),
  #                   stringsAsFactors = FALSE)
  #for(i in 1:length(x)){
  #  for(j in 1:length(x[[i]]@I)){
  #    for(k in 1:2){
  #      if(k == 1){
  #        data <- rbind(data,c(j, x[[i]]@I[j], "I", i))
  #      }
  #      else{
  #        data <- rbind(data,c(j, x[[i]]@R[j], "R", i))
  #      }
  #    }
  #  }
  #}
  #data$Population <- as.numeric(data$Population)
  #data$Index <- as.numeric(data$Index)
  #data$State <- as.factor(data$State)
  #data$ModelIndex <- as.factor(data$ModelIndex)
  #print(
  #  ggplot2::ggplot(data) +
  #    ggplot2::geom_line(
  #      ggplot2::aes(x = Index,
  #                   y = Population,
  #                   colour = State,
  #                   group = State:ModelIndex),
  #      alpha = 0.2) +
  #    ggplot2::scale_colour_manual(values = c("S" = "green",
  #                                            "I" = "red",
  #                                            "R" = "black"))
  #)
  #plot(summarise(x, average))
}
