
#COSTRUTTORI DI CLASSE per oggetti contenenti analisi
#' Checks if the object is an admitted 'sync analysis type'
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
is.sync = function(x){
  (inherits(x,"CCFBest")||inherits(x,"PMBest") )&& length(x)
}

CCFBest = function(sync, lag, ccf_matrix, lagSec, winSec, incSec, accelSec, weight_type, MAwin, MAinc,slopes)
{
  x = list("sync"=sync, "lag"=lag, "zero"=ccf_matrix["lag0"], "table"=ccf_matrix)
  class(x) = "CCFBest"
  attributes(x) = c(attributes(x), list(
      "lagSec"   = lagSec,
      "winSec"   = winSec,
      "incSec"   = incSec,
      "accelSec" = accelSec,
      "weight" = weight_type
  ))
  if(!missing(MAwin) && !missing(MAinc)){
    attributes(x) = c(attributes(x), list(
      "MAwinSec"   = MAwin,
      "MAincSec"   = MAinc
    ))
  }
  if(!missing(slopes)){
    attributes(x) = c(attributes(x), list(
      "slopes"   = TRUE
    ))
  }
  return(x)
}                    

#peakMatch
#' an AMICo analysis contains 3 objects:
#' a sync, and lag vector, containing synchrony and lag at the same frequency
#' of the original time series
#' and xBest, which is a table with the matched peaks positions and their similarity
#' also there are many attributes that save how the analysis was run
newAMICo = function(sync, lag, xBest, args)
{
  x = list("sync"=sync, "lag"=lag, "xBest"=xBest)
  class(x) = c("AMICo","PMBest") #@OBSOLETE pmbest
  attributes(x) = c(attributes(x), args)
  return(x) 
}


