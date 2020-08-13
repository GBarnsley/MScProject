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
  if(plot){
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
  if(plot){
    cat("\n",
        "The effective samples sizes are:","\n",
        " Beta: ")
    cat(NeffBeta)
    cat( "\n", " Gamma: ")
    cat(NeffGamma)
  }
  if(!plot){
    return(c(NeffBeta, NeffGamma))
  }
}
