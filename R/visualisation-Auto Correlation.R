#' A function that plots and prints ACF plots and effective sample size for samples from the
#' metropolis hastings algorithm using the acf() function.
#' @export
acfSamples <- function(samples, plot, timePoints){
  UseMethod("acfSamples", samples[[1]])
}
#' Method for SIR class, only plots for Beta and Gamma.
#' @export
acfSamples.SIR <- function(samples, plot = TRUE){
  Beta <- rep(NA, length(samples))
  Gamma <- rep(NA, length(samples))
  for(i in 1:length(samples)){
    Beta[i] <- samples[[i]]@Beta
    Gamma[i] <- samples[[i]]@Gamma
  }
  Beta <- acf(Beta, plot = FALSE)
  Gamma <- acf(Gamma, plot = FALSE)
  if(plot == TRUE){
    print(ggplot2::ggplot() +
            ggplot2::geom_segment(
              ggplot2::aes(
                x = Beta$lag,
                y = Beta$acf,
                yend = rep(0, length(Beta$lag)),
                xend = Beta$lag
            )) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = 0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = -0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::labs(
              x = "Lag",
              y = "ACF",
              title = "Beta"
            ))
    print(ggplot2::ggplot() +
            ggplot2::geom_segment(
              ggplot2::aes(
                x = Gamma$lag,
                y = Gamma$acf,
                yend = rep(0, length(Gamma$lag)),
                xend = Gamma$lag
            )) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = 0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = -0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::labs(
              x = "Lag",
              y = "ACF",
              title = "Gamma"
            ))
  }
  NeffBeta <- length(samples)/(1 + 2*sum(Beta$acf[-1]))
  NeffGamma <- length(samples)/(1 + 2*sum(Gamma$acf[-1]))
  cat("\n",
      "The effective samples sizes are:","\n",
      " Beta: ")
  cat(NeffBeta)
  cat( "\n", " Gamma: ")
  cat(NeffGamma)
}
#' Method for iSIR class, only plots for Beta, Gamma and newI.
#' @export
acfSamples.iSIR <- function(samples, plot = TRUE, timePoints = NA){
  Beta <- rep(NA, length(samples))
  Gamma <- rep(NA, length(samples))
  for(i in 1:length(samples)){
    Beta[i] <- samples[[i]]@Beta
    Gamma[i] <- samples[[i]]@Gamma
  }
  Beta <- acf(Beta, plot = FALSE)
  Gamma <- acf(Gamma, plot = FALSE)
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(samples)
  }
  newIs <- matrix(NA,
                  nrow = min(length(samples)-1, 10*log10(length(samples)) + 1),
                  ncol = length(timePoints))
  for(i in 1:ncol(newIs)){
    newI <- rep(NA, length(samples))
    for(j in 1:length(samples)){
      newI[j] <- samples[[j]]@newI[timePoints[i]]
    }
    newIs[,i] <- acf(newI, plot = FALSE)$acf
  }
  if(plot == TRUE){
    print(ggplot2::ggplot() +
            ggplot2::geom_segment(
              ggplot2::aes(
                x = Beta$lag,
                y = Beta$acf,
                yend = rep(0, length(Beta$lag)),
                xend = Beta$lag
              )) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = 0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = -0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::labs(
              x = "Lag",
              y = "ACF",
              title = "Beta"
            ))
    print(ggplot2::ggplot() +
            ggplot2::geom_segment(
              ggplot2::aes(
                x = Gamma$lag,
                y = Gamma$acf,
                yend = rep(0, length(Gamma$lag)),
                xend = Gamma$lag
              )) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = 0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = -0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::labs(
              x = "Lag",
              y = "ACF",
              title = "Gamma"
            ))
    for(i in 1:length(timePoints)){
      print(ggplot2::ggplot() +
              ggplot2::geom_segment(
                ggplot2::aes(
                  y = newIs[,i],
                  yend = rep(0, length(newIs[,i])),
                  x = 1:length(newIs[,i]),
                  xend = 1:length(newIs[,i])
                )
              ) +
              ggplot2::geom_hline(
                ggplot2::aes(
                  yintercept = 0.1
                ),
                alpha = 0.5,
                linetype = "dashed"
              ) +
              ggplot2::geom_hline(
                ggplot2::aes(
                  yintercept = -0.1
                ),
                alpha = 0.5,
                linetype = "dashed"
              ) +
              ggplot2::labs(
                x = "Lag",
                y = "ACF",
                title = paste(
                  "New I at time ",
                  timePoints[i]
                )
              ))
    }
  }
  NeffBeta <- length(samples)/(1 + 2*sum(Beta$acf[-1]))
  NeffGamma <- length(samples)/(1 + 2*sum(Gamma$acf[-1]))
  NeffTimePoints <- rep(NA, length(timePoints))
  for(i in 1:length(timePoints)){
    NeffTimePoints[i] <- length(samples)/(1 + 2*sum(newIs[-1,i]))
  }
  cat("\n",
      "The effective samples sizes are:","\n",
      " Beta: ")
  cat(NeffBeta)
  cat( "\n", " Gamma: ")
  cat(NeffGamma)
  for(i in 1:length(timePoints)){
    cat( "\n", " new I at time ")
    cat(timePoints[i])
    cat(": ")
    cat(NeffTimePoints[i])
  }
}
#' Method for rSIR class, only plots for Beta, Gamma and newI.
#' @export
acfSamples.rSIR <- function(samples, plot = TRUE, timePoints = NA){
  Beta <- rep(NA, length(samples))
  Gamma <- rep(NA, length(samples))
  for(i in 1:length(samples)){
    Beta[i] <- samples[[i]]@Beta
    Gamma[i] <- samples[[i]]@Gamma
  }
  Beta <- acf(Beta, plot = FALSE)
  Gamma <- acf(Gamma, plot = FALSE)
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(samples)
  }
  newRs <- matrix(NA,
                  nrow = min(length(samples)-1, 10*log10(length(samples)) + 1),
                  ncol = length(timePoints))
  for(i in 1:ncol(newRs)){
    newR <- rep(NA, length(samples))
    for(j in 1:length(samples)){
      newR[j] <- samples[[j]]@newR[timePoints[i]]
    }
    newRs[,i] <- acf(newR, plot = FALSE)$acf
  }
  if(plot == TRUE){
    print(ggplot2::ggplot() +
            ggplot2::geom_segment(
              ggplot2::aes(
                x = Beta$lag,
                y = Beta$acf,
                yend = rep(0, length(Beta$lag)),
                xend = Beta$lag
              )) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = 0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = -0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::labs(
              x = "Lag",
              y = "ACF",
              title = "Beta"
            ))
    print(ggplot2::ggplot() +
            ggplot2::geom_segment(
              ggplot2::aes(
                x = Gamma$lag,
                y = Gamma$acf,
                yend = rep(0, length(Gamma$lag)),
                xend = Gamma$lag
              )) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = 0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::geom_hline(
              ggplot2::aes(
                yintercept = -0.1
              ),
              alpha = 0.5,
              linetype = "dashed"
            ) +
            ggplot2::labs(
              x = "Lag",
              y = "ACF",
              title = "Gamma"
            ))
    for(i in 1:length(timePoints)){
      print(ggplot2::ggplot() +
              ggplot2::geom_segment(
                ggplot2::aes(
                  y = newRs[,i],
                  yend = rep(0, length(newRs[,i])),
                  x = 1:length(newRs[,i]),
                  xend = 1:length(newRs[,i])
                )
              ) +
              ggplot2::geom_hline(
                ggplot2::aes(
                  yintercept = 0.1
                ),
                alpha = 0.5,
                linetype = "dashed"
              ) +
              ggplot2::geom_hline(
                ggplot2::aes(
                  yintercept = -0.1
                ),
                alpha = 0.5,
                linetype = "dashed"
              ) +
              ggplot2::labs(
                x = "Lag",
                y = "ACF",
                title = paste(
                  "New R at time ",
                  timePoints[i]
                )
              ))
    }
  }
  NeffBeta <- length(samples)/(1 + 2*sum(Beta$acf[-1]))
  NeffGamma <- length(samples)/(1 + 2*sum(Gamma$acf[-1]))
  NeffTimePoints <- rep(NA, length(timePoints))
  for(i in 1:length(timePoints)){
    NeffTimePoints[i] <- length(samples)/(1 + 2*sum(newRs[-1,i]))
  }
  cat("\n",
      "The effective samples sizes are:","\n",
      " Beta: ")
  cat(NeffBeta)
  cat( "\n", " Gamma: ")
  cat(NeffGamma)
  for(i in 1:length(timePoints)){
    cat( "\n", " new R at time ")
    cat(timePoints[i])
    cat(": ")
    cat(NeffTimePoints[i])
  }
}
