##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  
## DyadFilter.R
# FILTERS!, and other time series processing goodness
#
############################################################################################################ 
## changelog
# v1.1 - added setArtefacts function to remove bad parts of signal
# v1.0 - fully 'stream 1.2' compliant. All ts ARE ok, with explicit start calls.
#         
#
############################################################################################################ 
## ToDo
# -all filters should apply on DyadStreams only. The Experiment/signal stuff 
#  should be managed by wrapper functions
#
############################################################################################################ 

#' setArtefacts
#' setArtefacts legge una tabella o lista con le componenti chiamate esattamente: 'session', 'start' ed 'end'
#' e mette ad NA tali sequenze nel segnale specificato.
#' end può anche essere una stringa che contiene le parole 'fine' o 'end'
#' 
#' @param x 
#' @param startEnd 
#' @param signal 
#'
#' @return
#' @export
#'
#' @examples
setArtefacts <- function(x, startEnd, signal="SC"){
  names(startEnd) = tolower(names(startEnd))
  if("factor" %in% sapply(startEnd,class)) stop ("factors not supported in startEnd")
  UseMethod("setArtefacts",x)
}

#' @export
setArtefacts.DyadExperiment <- function(x, startEnd, signal="SC") {
  if(is.data.frame(startEnd)){
    if(ncol(startEnd)!=3 || !all.equal(names(startEnd), c('session', 'start', 'end') ))
      stop("startEnd must have 3 columns: 'session', 'start', 'end' of equal size")
    sel = startEnd
  } else if(is.list(startEnd)) {
    if(length(startEnd)!=3 || var(sapply(startEnd, length))!=0 )
      stop("startEnd must have 3 vectors 'session', 'start', 'end' of equal size")
    sel = as.data.frame(startEnd,stringsAsFactors = F)
  } else stop("startEnd must be a list or dataframe")
  
  for(i in unique(sel$session) ){
    cat("\r\ncleaning session:",i,"\r\n")
    miniSel = sel[sel$session == i, 2:3]
    x[[i]][[signal]] = setArtefacts(x[[i]][[signal]],miniSel,signal)
  }
  x
}

#' @export
setArtefacts.DyadSignal <- function(x, startEnd, signal="SC") {
  #1 controlla validità di startEnd
  if(!is.data.frame(startEnd)){
    if(length(startEnd)!=2 || length(startEnd[[1]])!=length(startEnd[[2]]))
      stop("startEnd must have 2 vectors 'start', 'end' of equal size")
    sel = as.data.frame(startEnd,stringsAsFactors = F)
  } else sel = startEnd
  if(!all.equal(names(sel),c('start', 'end')) ) stop("startEnd names must be 'start','end'")
  ref = 0
  duration = duration(x)
  print(sel)
  
  for(i in 1:nrow(sel)){
    a = timeMaster(sel[i,1],out="sec") -  start(x)[1]
    if( grepl("end|fine",sel[i,2],ignore.case = T)  ) b = duration
    else b = timeMaster(sel[i,2],out="sec") -  start(x)[1]
    if(a>duration || b > duration) stop ("one start or end were bigger than the signal length")
    if(b<a) stop ("'start' cannot be greater than 'end' ")
    if(a < ref || b < ref) stop("all start and end definition must be progressively increasing and greater than signal start value")
    ref = b
    sel[i,1] = a * sampRate(x)
    if(b == duration) sel[i,2] =  length(x$valid) #elminina tutto fino alla fine
    else  sel[i,2] = b * sampRate(x)
  }

  #3 sostituisci con NA i segmenti
  for(i in 1:nrow(sel)){
    x$valid[sel[i,1]:sel[i,2]] = FALSE
  }
  
  #4 aggiorna i metadati
  attributes(x)["filter"] = paste0(attr(x,"filter"), "--> NAartifact")
  x
}



############################################################################################################
############################################################################################################
############################################################################################################
## Decimation and interpolation functions

#' signalDecimate
#' Takes an objects of class DyadSignal and downsamples it by taking a sample each oldSampRate / newSampRate
#'
#' @param signal 
#' @param newSampRate 
#'
#' @return
#' @export
#'
#' @examples
signalDecimate = function (signal, newSampRate) {
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  if (sampRate(signal) <= newSampRate) 
    stop("decimate only downsamples! newSampRate:",newSampRate," old:",sampRate(signal),call. = F)
  if (sampRate(signal) %% newSampRate != 0) 
    stop("newSampRate must be an integer divisor of old sampRate ! newSampRate:",newSampRate," old:",sampRate,call. = F)
  q = floor(sampRate(signal) / newSampRate)  
  resPat = ts(signal$s1[seq(1,   length(signal$s1)  , by = q)], start=start(signal$s1),   frequency=newSampRate)
  resCli = ts(signal$s2[seq(1, length(signal$s2), by = q)], start=start(signal$s2), frequency=newSampRate)
  resVal = ts(signal$valid[seq(1, length(signal$valid), by = q)], start=start(signal$valid), frequency=newSampRate)
  
  signal$s1   = cloneDyadStream(resPat, signal$s1)
  signal$s2 = cloneDyadStream(resCli, signal$s2)
  signal$valid = resVal
  signal$time = time(signal$s1)
  attr(signal,"sampRate") = newSampRate
  return(signal)
}

### onlyy used one in connect area of graphic library. maybe delete, maybe integrate
# streamDecimate = function (x, newSampRate) {
#   if(!is(x,"DyadStream")) stop("Only objects of class DyadStream can be processed by this function")
#   if (frequency(x) <= newSampRate) 
#     stop("decimate only downsamples! newSampRate:",newSampRate," old:",frequency(x),call. = F)
#   if (frequency(x) %% newSampRate != 0) 
#     stop("newSampRate must be an integer divisor of old sampRate ! newSampRate:",newSampRate," old:",frequency(x),call. = F)
#   q = floor(frequency(x) / newSampRate)  
#   res = ts(x[seq(1, length(x), by = q)], start=start(x), frequency=newSampRate)
#   cloneDyadStream(res, x)
# }

#Interpolate windowed data to original samplerate. 
#useful to overlay computed indexes on original timeseries
winInter = function(windowsList, winSec, incSec, sampRate){
  if(class(windowsList)!="list") {windowsList = list(windowsList)}
  #cat("Interpolating",incSec*sampRate,"sample between each HRV datapoint (linear) \r\n")
  inc=incSec*sampRate
  win= winSec*sampRate
  nList=length(windowsList)
  Map(function(daba,i){
    #prog(i,nList)
    data.frame(apply(daba,2,function (series){
      approx(x    = seq(1,(length(series)*inc), by=inc)+ceiling(win/2)-1, #windowed datapoints must be at the center of the window
             y    = series, #the exact calculated values
             xout = seq(1,(length(series)-1)*inc + win,by=1) #all samples covered by the windows
      )$y
    }))
  },windowsList,seq_along(windowsList))
}



############################################################################################################
############################################################################################################
############################################################################################################
#### Moving average filters

#movAv has improved managing of burn-in and burn-out sequences
#if in doubt use this.
#' Title
#'
#' @param x 
#' @param win if x has a frequency, it will be the window size in seconds. Otherwise in samples.
#' @param remove 
#'
#' @return
#' @details please note that MA filters are not adequate to smooth sync calculation that shifts
#' often between extreme negative and positive values.
#' @export
#'
#' @examples
#' a = c(rep(0,20),rep(5,20),rep(10,20),rep(20,20))
#' plot(a)
#' lines(movAv(a,10,F))
#' lines(movAv(a,20,F),col=3)
#' b=DyadStream(a,"test",col=2,frequency = 10, start=0)
#' lines(movAv(b,2,F),col=color(b),lty=3)
#' sig=list("s1" = b, s2 = b)
#' class(sig) ="DyadSignal" 
#' lines(movAv(sig$s1,2,F),col="skyblue",lwd=2)


movAv = function(x, win, remove=F) {
  UseMethod("movAv",x)
}
#' @export
movAv.DyadSignal = function(signal,win,remove=F){
  res = lapply(list(signal$s1,signal$s2), movAv, win=win, remove=remove)
  signal$s1 = res[[1]]
  signal$s2 = res[[2]]
  return(signal)
}
#' @export
movAv.DyadStream = movAv.ts = function(x,win,remove=F){
  sampRate=frequency(x)
  y = movAv.default(x,win*sampRate,remove=F)
  classAttr(y) = classAttr(x)
  y
}
#' @export
movAv.default = function(x,win,sampRate,remove=F){
  win2 = floor(win/2)
    len = length(x)
    f = unlist(lapply(seq_along(x),function(t){
      i1 = ifelse(t <= win2, 1, (t-win2) )
      i2 = ifelse(t+win2 >= len, len, (t+win2) )
      sum(x[i1:i2])/length(i1:i2)
    }))
    if(remove){
      x-f
    }else  f
}


## znorm normalizes a time series. Good to compare with other TS
#' Title
#'
#' @param a a DyadStream object or a numeric vector to be passed to ts()
#'
#' @return a DyadStream or a ts object
#' @export
#'
#' @examples
znorm = function(a){
  if(is.DyadStream(a))
    cloneDyadStream(ts(scale(a)[,1],frequency=frequency(a),start=start(a),end=end(a)),a)
  else 
    ts(scale(a)[,1],frequency=frequency(a),start=start(a))
}

## this function removes the moving average of a signal, but using nonoverlapping windows.
## this is useful to plot in a horizontal panel a signal with huge trends.
stepCenter = function(a, winSec=60){
  if(is.DyadStream(a)) x=a
  
  freq = frequency(a)
  winSam = winSec*freq
  n_win = ceiling((length(a)-winSam+1)/winSam)
  resid = length(a)-n_win*winSam
  for(i in 0:(n_win-1)){
    al = (i*winSam +1)
    bl = (i*winSam +winSam)
    a[al:bl] = a[al:bl] - mean(a[al:bl])
  }
  if (length(resid)>0)
    a[(length(a)-resid+1):length(a)] = a[(length(a)-resid+1):length(a)] - mean(a[(length(a)-resid+1):length(a)])
  ifelse(is.DyadStream(x), cloneDyadStream(a,x) , a )
}



