
#COSTRUTTORI DI CLASSE per oggetti contenenti analisi
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
PMBest = function(sync, lag, xBest, lagSec, sgol_p, sgol_n, weightMalus)
{
  x = list("sync"=sync, "lag"=lag, "xBest"=xBest)
  class(x) = "PMBest"
  attributes(x) = c(attributes(x), list(
      "lagSec"   = lagSec,
      "sgol_p"   = sgol_p,
      "sgol_n"   = sgol_n,
      "weightMalus" = weightMalus
  ))
  return(x) 
}
