#     ____  _               _                                   
#    / __ \| |             | |                                  
#   | |  | | |__  ___  ___ | | ___  ___  ___ ___ _ __   ___ ___ 
#   | |  | | '_ \/ __|/ _ \| |/ _ \/ __|/ __/ _ \ '_ \ / __/ _ \
#   | |__| | |_) \__ \ (_) | |  __/\__ \ (_|  __/ | | | (_|  __/
#    \____/|_.__/|___/\___/|_|\___||___/\___\___|_| |_|\___\___|
#                                                               
#                                                               
#' Here goes deprecated and obsolete functions calls, that will eventually get deleted

#' @rdname extractEpochs
#' @export
catExtractLong = function(x, signal="", epochStream="",by="",FUN=""){stop("this function has been renamed to extractEpochs ")}



## legacyPeakFinder is an old version of peakFinder not implementing max time 
legacyPeakFinder = function(x, sgol_p = 2, sgol_n = 25, mode=c("peaks","valleys","both"), correctionRangeSeconds = 0.5){
  sampRate = frequency(x)
  smooth_x = signal::sgolayfilt(x,  p =sgol_p, n = sgol_n, m=0)
  fdderiv1  = diff(smooth_x)
  mode = match.arg(mode)
  pik = sign(embed(fdderiv1,2)) #embed appaia al segnale il segnale laggato di 1 samp
  s = pik[,2] - pik[,1] #that's where the magic happens
  
  if(mode=="peaks") {
    pikboo = c(s ==  2, FALSE) #embed perde 1 sample, quindi aggiungi un FALSE alla fine
    piksam = which(pikboo) #in quali sample c'è un'inversione di segno della derivata?
    pv = rep("p",length(piksam))
  } else if(mode=="valleys") {
    pikboo = c(s == -2, FALSE) 
    piksam = which(pikboo) 
    pv = rep("v",length(piksam))
  } else if(mode=="both") {
    pikboo = c(abs(s) ==  2, FALSE)
    piksam = which(pikboo) 
    s = s[s!=0]
    pv = s
    pv[s>0] = "p"
    pv[s<0] = "v"
  }
  
  #correzione manuale: cerca il valore più alto nei dintorni del cambio di derivata
  for(v in seq_along(piksam)){
    #individua il range con 0.5s prima e dopo la valle della derivata (esclusi gli estremi inzio e fine ts)
    search_interval = x[
      max(1,piksam[v]-correctionRangeSeconds*sampRate):min(length(x),piksam[v]+correctionRangeSeconds*sampRate)
    ]
    #trova il più piccolo e aggiorna pikmsam
    if(mode=="peaks")
      piksam[v] = max(1, piksam[v]+  (which.max(search_interval) - round(length(search_interval)/2)))
    else if (mode=="valleys")
      piksam[v] = max(1, piksam[v]+  (which.min(search_interval) - round(length(search_interval)/2)))
    else
      piksam[v] = max(1, piksam[v]+  (which.minmax(search_interval,pv[v]) - round(length(search_interval)/2)))
    
    piksam[v] = min(piksam[v], length(x))
  }
  #ricrea pikmboo dai sample corretti
  pikboo = rep(F,length(pikboo))
  pikboo[piksam] = T
  piks = time(x)[piksam]
  list("bool" = pikboo,
       "samples" = piksam,
       "seconds" = piks,
       "type" = pv)
}


#vecchio "get"
#' Title
#'
#' @param x 


#' @export
start.DyadSignal <- function (x) {attr(x,"start")}
#' @export
end.DyadSignal <- function (x) {attr(x,"end")}
#' @export

#' @export
s1Name <- function(x) {attr(x,"s1Name")}
#' @export
s2Name <- function(x) {attr(x,"s2Name")}

#' @export
sessionId = function(dyadSession){
  attr(dyadSession,"sessionId")
}
#' @rdname sessionId
#' @export

#' @export
dyadId = function(dyadSession){
  attr(dyadSession,"dyadId")
} 
#' @export


## ## GROUP ID
## get
#' @export
groupId = function(x) {
  UseMethod("groupId", x)
}

#' @export
groupId.DyadExperiment = function(x){
  `Group Ids` = sapply(x, attr, "groupId")
  table(`Group Ids`)
}

#' @export
groupId.DyadSession = function(x){
  attr(x,"groupId")
}

#' @export
name <- function(x) {
  attr(x,"name")
}

## set
#' @export
`groupId<-` = function(x, value){
  UseMethod("groupId<-",x)
}
#' @export
`groupId<-.DyadExperiment` <- function(x, value){
  namesx = character(length(x))
  for(i in 1:length(x)){
    groupId(x[[i]]) = value
    namesx[[i]] = name(x[[i]])
  }
  names(x) = namesx
  x
}
#' @export
`groupId<-.DyadSession` <- function(x, value){
  name = paste(value,lead0(sessionId(x)),dyadId(x),sep="_")
  attr(x, "groupId") <- value
  attr(x, "name")    <- name
  x
}
#' 
#' 
#' 
#' #' @export
#' uniqueId.DyadExperiment = function(x){
#'   `Group Ids` = sapply(x, attr, "groupId")
#'   table(`Group Ids`)
#' }
#' 
#' #' @export
#' uniqueId.DyadSession = function(x){
#'   attr(x,"groupId")
#' }
#' 
#' 
#' 
#' 
#' 
#' 

#' @export
sampRate <-function(x) {
  if(is.ts(x)) frequency(x)
  else attr(x,"sampRate")
}

#' @export
frequency.DyadSignal = function(x){attr(x,"sampRate")}
