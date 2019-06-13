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


#' apply a function on each dyad's member signal
#'
#' @param signal DyadSignal Object 
#' @param FUN 
#' @param ... additional parameters of FUN
#' @param newAttributes named list of attributes to be added or replaced
#'
#' @return
#' @export
#'
#' @examples
#' 
signalFilter = function (x, FUN, ..., newAttributes=NULL, signals="all") {
  UseMethod("signalFilter",x)
}
#' @export
signalFilter.DyadExperiment = function (x, FUN, ..., newAttributes=NULL, signals="all") {
  #for each session
  experiment2 = Map(function (session,nSession){
    #find signal names
    if(length(signals)==1 && signals=="all") {
      sigs = names(session)[sapply(session,is.DyadSignal)]
    } else sigs = signals
    #apply signalFilter on each DyadSignal object
    session[names(session) %in% sigs] = lapply(session[names(session) %in% sigs],  signalFilter, FUN, ..., newAttributes, signals=NULL)
    prog(nSession,length(x))
    return(session)
  },x,seq_along(x))
  experiment2 = cloneAttr(x,experiment2)
  return(experiment2)
}
#' @export
signalFilter.DyadSignal = function (x, FUN, ..., newAttributes=NULL, signals=NULL) {
  if(!is(x,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  FUN = match.fun(FUN)
  ress1 = FUN(x$s1, ...)
  ress2 = FUN(x$s2, ...)
  if(!is.ts(ress1)) ress1 = ts(ress1, start=start(x$s1), frequency=frequency(x$s1))
  if(!is.ts(ress2)) ress2 = ts(ress2, start=start(x$s2), frequency=frequency(x$s2))
  
  x$s1 = cloneAttr(x$s1, ress1)
  x$s2 = cloneAttr(x$s2, ress2)
  x$time = time(x$s1)
  if(!is.null(newAttributes)){
    if("filter"%in%names(newAttributes)){
      newAttributes[["filter"]] = paste0(attributes(x)[["filter"]]," -> ",newAttributes[["filter"]])
    }
    attributes(x)[names(newAttributes)] = newAttributes
  }
  return(x)
}


#' setArtefacts
#' setArtefacts legge una tabella o lista con le componenti chiamate esattamente: 'dyad', 'session', 'start' ed 'end'
#' salva le informazioni corrispondenti in una tabella $artefacts nell'oggetto DyadSignal e mette a FALSE le epoche
#' corrispondenti nello stream 'valid' dell'oggetto.
#' 'end' può anche essere una stringa che contiene le parole 'fine' o 'end'
#' 
#' @param x a Dyad... object
#' @param startEnd a data.frame (or list) with the following components: 'dyad', 'session', 'start', 'end'.
#' @param signal string specifying the name of the signal
#'
#' @return
#' @export
#'
#' @examples
setArtefacts <- function(x, startEnd, signal){
  names(startEnd) = tolower(names(startEnd))
  if("factor" %in% sapply(startEnd,class)) stop ("factors not supported in startEnd")
  UseMethod("setArtefacts",x)
}

#' @export
setArtefacts.DyadExperiment <- function(x, startEnd, signal) {
  if(is.data.frame(startEnd)){
    if(ncol(startEnd)!=4 || !all.equal(names(startEnd), c('dyad','session', 'start', 'end') ))
      stop("startEnd must have 4 columns: 'dyad', 'session', 'start', 'end' of equal size")
    sel = startEnd
  } else if(is.list(startEnd)) {
    if(length(startEnd)!=4 || var(sapply(startEnd, length))!=0 )
      stop("startEnd must have 3 vectors 'session', 'start', 'end' of equal size")
    sel = as.data.frame(startEnd,stringsAsFactors = F)
  } else stop("startEnd must be a list or dataframe")
  
  for(j in unique(sel$dyad) ){
    for(i in unique(sel$session) ){
      listKey = which(sapply(x,sessionId)==i & sapply(x,dyadId)==j) #questo è importante per selezionare la seduta giusta
      if(length(listKey)>1) stop("multiple matches found in session:",j,lead0(i))
      if(length(listKey)==1){
        cat("\r\ncleaning session:",j,lead0(i),"\r\n")
        miniSel = sel[sel$dyad == j & sel$session == i, 3:4]
        x[[listKey]][[signal]] = setArtefacts(x[[listKey]][[signal]],miniSel,signal)
      }
      }}
  x
}

#' @param x 
#'
#' @param startEnd 
#' @param signal 
#'
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
  for(i in 1:nrow(sel)){
    sel[grepl("end|fine",sel[,2],ignore.case = T),2] = timeMaster(duration,"min")
  }
  # cat("\r\n ! - ",str(sel))
  
  # if(nrow(sel)==1)
  x$artefacts = data.frame(numeric(nrow(sel)))
  x$artefacts$start = timeMaster(sel[,1], out="s")
  x$artefacts$end   = timeMaster(sel[,2], out="s")
  x$artefacts[1] = NULL
  # print(str(x$artefacts))

  #questa roba è qui per compatibilità. Idealmente usa solo la tabella artefacts ##############
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
  for(i in 1:nrow(sel)){  #sostituisci con NA i segmenti
    x$valid[sel[i,1]:sel[i,2]] = FALSE
  }
  ###############################
  
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
  ress1 = ts(signal$s1[seq(1,   length(signal$s1)  , by = q)], start=start(signal$s1),   frequency=newSampRate)
  ress2 = ts(signal$s2[seq(1, length(signal$s2), by = q)], start=start(signal$s2), frequency=newSampRate)
  resVal = ts(signal$valid[seq(1, length(signal$valid), by = q)], start=start(signal$valid), frequency=newSampRate)
  
  signal$s1 = cloneAttr(signal$s1, ress1)
  signal$s2 = cloneAttr(signal$s2, ress2)
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
#' @param win 
#' @param inc 
#' @param sampRate 
#' @param remove 
#'
#' @return
#' @export
#' 
#' 
movAv = function(x,win,inc,sampRate=frequency(x)){
  warning("this function may be broken. inc is not used")
  n = winSec*sampRate
  b = unlist(lapply(seq_along(a), function(t){
    i1 = t-(n-1)/2; if(i1<1) i1=1;
    i2 = t+(n-1)/2;
    res = sum(a[i1:i2])/n
    res[is.na(res)]=0
    return(res)
  }))
  ts(b,start=start(a)[1]+winSec/2,frequency = sampRate)
  
}
# movAv = function(x,win,inc,sampRate=frequency(x)){
#   stop("this function is broken")
#   win = win*sampRate
#   inc = inc*sampRate
#   win2 = floor(win/2)
#   len = length(x)
#   
#   
#   n_win = ceiling((length(x)-win+1)/inc)
#   res = numeric(n_win)
#  
#   for(i in 1:n_win-1){
#     a = (i*inc +1)
#     b= (i*inc +win)
#     res[i] = mean(x[a:b],na.rm=T)
#     
#   }
#   ts(res,  start=c(floor(time(x)[win2]),cycle(x)[win2]), frequency = inc)
# 
# }

movAvSlope = function(x,win,inc,sampRate=frequency(x)){
  warning("this function may be broken")
  win = win*sampRate
  inc = inc*sampRate
  n_win = ceiling((length(x)-win+1)/inc)
  res = numeric(n_win)
  for(i in 1:nwin-1){
    a = (i*inc +1)
    b= (i*inc +win)
    res[i] = (x[b]-x[a])/win
  }
  ts(res, frequency = inc, start=start(x))
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
    cloneAttr(a,ts(scale(a)[,1],frequency=frequency(a),start=start(a),end=end(a)))
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



