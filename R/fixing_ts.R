#' ts method for scale()
#' @export

scale.ts = function(x, center = TRUE, scale = TRUE){
  asd = attributes(x)
  x = NextMethod(x)
  attributes(x) = c(attributes(x),asd)
  return(x)
}

#' Time Series Start
#' this returns the start of a ts object in time units
#'
#' @param x a ts object
#' @export
tss = function(x){tsp(x)[1]}
#' Time Series End
#' this returns the end of a ts object in time units
#'
#' @param x a ts object
#' @export
tse = function(x){tsp(x)[3]}