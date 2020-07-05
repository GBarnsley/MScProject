#' Probability Generator
#'
#' Converts a given value into a value in the region (0,1).
#' Uses the exponential cdf. This funtion is only defined
#' seperately so that it can be easily changed if needed.
#'
#' @param x some non-negative value
#' @return a value on the range (0,1)
#' @export
probGen <- nimble::nimbleFunction(
  run = function(x = double(0)) {
    returnType(double(0))
    return(1 - exp(-x))
  }
)
