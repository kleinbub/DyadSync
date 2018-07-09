SYNC_CLASSES = c("CCFBest","PMBest")

#' hasSync
#' Recursive function which checks if any of the nested object(s) (typically a DyadSignal)
#' contains at least one children of a class defined in the global constant SYNC_CLASSES
#' @param x 
#' @export
hasSync = function(x){
  res = unique(rapply(x,class))
  if(any(res%in%SYNC_CLASSES)){
    which(res%in%SYNC_CLASSES)
  } else FALSE
}

# JUST WHY?!?
# #' settings
# #'
# #' @param x 
# #' @param which 
# #'
# #' @return
# #' @export

# settings = function (x, which)
# {
#   if(missing(which)){
#     attr(x,"settings")
#   } else attr(x,"settings")[[which]]
# }



#COSTRUTTORI DI CLASSE per oggetti contenenti analisi

CCFBest = function(sync, lag, ccf_matrix, lagSec, winSec, incSec, accelSec, weight)
{
  x = list("sync"=sync, "lag"=lag, "table"=ccf_matrix)
  class(x) = "CCFBest"
  attributes(x) = c(attributes(x), list(
      "lagSec"   = lagSec,
      "winSec"   = winSec,
      "incSec"   = incSec,
      "accelSec" = accelSec,
      "weight" = weight
  ))
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
