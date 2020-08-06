#' A function that plots and prints ACF plots and effective sample size for samples from the
#' metropolis hastings algorithm using the acf() function.
#' @export
acfSamples <- function(epiModel, plot){
  UseMethod("acfSamples")
}
#' Method for SIR class, only plots for Beta and Gamma.
#' @export
acfSamples.SIR <- function(epiModel, plot = TRUE){
  Beta <- acf(epiModel@Samples[,'Beta'], plot = FALSE)
  Gamma <- acf(epiModel@Samples[,'Gamma'], plot = FALSE)
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
  NeffBeta <- nrow(epiModel@Samples)/(1 + 2*sum(Beta$acf[-1]))
  NeffGamma <- nrow(epiModel@Samples)/(1 + 2*sum(Gamma$acf[-1]))
  cat("\n",
      "The effective samples sizes are:","\n",
      " Beta: ")
  cat(NeffBeta)
  cat( "\n", " Gamma: ")
  cat(NeffGamma)
}
#' Method for iSIR class, only plots for Beta, Gamma and newI.
#' @export
acfSamples.iSIR <- function(epiModel, plot = TRUE, timePoints = NA){
  Beta <- acf(epiModel@Samples[,'Beta'], plot = FALSE)
  Gamma <- acf(epiModel@Samples[,'Gamma'], plot = FALSE)
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(epiModel)
  }
  newIs <- list()
  for(i in 1:length(timePoints)){
    newIs[[i]] <- acf(epiModel@Samples[,paste0("newI[",timePoints[i],"]")], plot = FALSE)$acf
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
                  y = newIs[[i]],
                  yend = rep(0, length(newIs[[i]])),
                  x = 1:length(newIs[[i]]),
                  xend = 1:length(newIs[[i]])
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
  NeffBeta <- nrow(epiModel@Samples)/(1 + 2*sum(Beta$acf[-1]))
  NeffGamma <- nrow(epiModel@Samples)/(1 + 2*sum(Gamma$acf[-1]))
  NeffTimePoints <- rep(NA, length(timePoints))
  for(i in 1:length(timePoints)){
    NeffTimePoints[i] <- nrow(epiModel@Samples)/(1 + 2*sum(newIs[[i]][-1]))
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
acfSamples.rSIR <- function(epiModel, plot = TRUE, timePoints = NA){
  Beta <- acf(epiModel@Samples[,'Beta'], plot = FALSE)
  Gamma <- acf(epiModel@Samples[,'Gamma'], plot = FALSE)
  if(identical(timePoints, NA)){
    timePoints <- generateTimepoints(epiModel)
  }
  newIs <- list()
  for(i in 1:length(timePoints)){
    newIs[[i]] <- acf(epiModel@Samples[,paste0("newR[",timePoints[i],"]")], plot = FALSE)$acf
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
                  y = newIs[[i]],
                  yend = rep(0, length(newIs[[i]])),
                  x = 1:length(newIs[[i]]),
                  xend = 1:length(newIs[[i]])
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
  NeffBeta <- nrow(epiModel@Samples)/(1 + 2*sum(Beta$acf[-1]))
  NeffGamma <- nrow(epiModel@Samples)/(1 + 2*sum(Gamma$acf[-1]))
  NeffTimePoints <- rep(NA, length(timePoints))
  for(i in 1:length(timePoints)){
    NeffTimePoints[i] <- nrow(epiModel@Samples)/(1 + 2*sum(newIs[[i]][-1]))
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
