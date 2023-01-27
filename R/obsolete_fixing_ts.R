#' #' ts method for scale()
#' #' @export
#' 
#' scale.ts = function(x, center = TRUE, scale = TRUE){
#'   asd = attributes(x)
#'   x = NextMethod(x)
#'   attributes(x) = c(attributes(x),asd)
#'   return(x)
#' }
#' 
#' #' Time Series Start
#' #' this returns the start of a ts object in time units
#' #'
#' #' @param x a ts object
#' #' @export
#' tss = function(x){UseMethod("tss",x)}
#' #' @export
#' tss.ts = function(x){return(tsp(x)[1L])}
#' #' @export
#' tss.DyadSignal = function(x){return(start(x)[1L])}
#' 
#' #start(x)[1]*frequency(x)+start(x)[2]-1
#' 
#' #' Time Series End
#' #' this returns the end of a ts object in time units
#' #'
#' #' @param x a ts object
#' #' @export
#' tse = function(x){UseMethod("tse",x)}
#' #' @export
#' tse.ts = function(x){return(tsp(x)[2L])}
#' # end(x)[1]*frequency(x)+end(x)[2]-1
#' #' @export
#' tse.DyadSignal = function(x){return(end(x)[1L])}
#' 
#' #' Duration of a ts object in time units
#' #' @param x a ts object
#' #'
#' #' @export
#' duration = function(x){UseMethod("duration",x)}
#' #' @export
#' duration.ts = function(x){return(tse(x)-tss(x))}
#' 
#' 
#' 
