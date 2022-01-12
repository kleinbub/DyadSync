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


#' applies a function on the s1 and s2 objects of a DyadSignal object
#'
#' @param x a DyadExperiment, DyadSession, or DyadSignal Object 
#' @param FUN 
#' @param newAttributes named list of attributes to be set on the resulting object. "filter" attributes will be added instead.
#' @param signals a vector of strings defining the names of the signals on which to run the filter
#' @param ... additional arguments of FUN
#' @return a DyadExperiment, DyadSession, or DyadSignal Object with the filtered signals
#' @export
#'

signalFilter = function (x, FUN, newAttributes=NULL, signals="all", ...) {
  UseMethod("signalFilter",x)
}
#' @export
signalFilter.DyadExperiment = function (x, FUN, newAttributes=NULL, signals="all", ...) {
  #for each session
  experiment2 = Map(function (session,nSession){
    session = signalFilter(session, FUN=FUN, newAttributes=newAttributes, signals=signals, ...)
    prog(nSession,length(x))
    return(session)
  },x,seq_along(x))
  experiment2 = cloneAttr(x,experiment2)
  return(experiment2)
}

#' @export
signalFilter.DyadSession = function (x, FUN, newAttributes=NULL, signals="all", ...) {
  #find signal names
  if(length(signals)==1 && signals=="all") {
    sigs = names(x)[sapply(x,is.DyadSignal)]
  } else sigs = signals
  #apply signalFilter on each DyadSignal object
  x[names(x) %in% sigs] = lapply(x[names(x) %in% sigs], function(x){signalFilter(x, FUN=FUN, newAttributes=newAttributes, ...)})
  x
}

#' @export
signalFilter.DyadSignal = function (x, FUN, newAttributes=NULL, signals=NULL, ...) {
  FUN = match.fun(FUN)
  ress1 = FUN(x$s1, ...)
  ress2 = FUN(x$s2, ...)
  if(!is.ts(ress1)) ress1 = ts(ress1, start=start(x$s1), frequency=frequency(x$s1))
  if(!is.ts(ress2)) ress2 = ts(ress2, start=start(x$s2), frequency=frequency(x$s2))
  
  x$s1 = cloneAttr(x$s1, ress1)
  x$s2 = cloneAttr(x$s2, ress2)
  # x$time = time(x$s1)
  
  attr(x,"start") = start(x$s1)
  attr(x,"end")   = end(x$s1)
  attr(x,"sampRate") = frequency(x$s1)

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
  #1 controlla validità di startEnd
  names(startEnd) = tolower(names(startEnd))
  if("factor" %in% sapply(startEnd,class)) stop ("factors not supported in startEnd")
  if(is.data.frame(startEnd)){
    if(ncol(startEnd)!=4 || !all.equal(names(startEnd), c('dyad','session', 'start', 'end') ))
      stop("startEnd must have 4 columns: 'dyad', 'session', 'start', 'end' of equal length")
    sel = startEnd
  } else if(is.list(startEnd)) {
    if(length(startEnd)!=4 || var(sapply(startEnd, length))!=0 )
      stop("startEnd must have 3 vectors 'session', 'start', 'end' of equal length")
    sel = as.data.frame(startEnd,stringsAsFactors = F)
  } else stop("startEnd must be a list or dataframe")

  nClean = 0  
  for(j in unique(sel$dyad) ){
    for(i in unique(sel$session) ){
      listKey = which(sapply(x,sessionId)==i & sapply(x,dyadId)==j) #questo è importante per selezionare la seduta giusta
      if(length(listKey)>1) stop("multiple matches found in session:",j,lead0(i))
      if(length(listKey)==1){
        cat("\r\ncleaning session:",j,lead0(i),"\r\n")
        miniSel = sel[sel$dyad == j & sel$session == i, ]
        x[[listKey]][[signal]] = setArtefacts(x[[listKey]][[signal]],miniSel)
        nClean = nClean +1
      }
    }}
  if(nClean == 0) warning("no sessions were affected. Maybe check dyadId and sessionId in both DyadExperiment object and startEnd dataset")
  x
}

#' @param x 
#'
#' @param startEnd 
#'
#' @export
setArtefacts.DyadSignal <- function(x, startEnd) {
  
  #1 controlla validità di startEnd
  names(startEnd) = tolower(names(startEnd))
  if("factor" %in% sapply(startEnd,class)) stop ("factors not supported in startEnd")
  if(is.data.frame(startEnd)){
    if(ncol(startEnd)!=4 || !all.equal(names(startEnd), c('dyad','session', 'start', 'end') ))
      stop("startEnd must have 4 columns: 'dyad', 'session', 'start', 'end' of equal length")
    sel = startEnd
  } else if(is.list(startEnd)) {
    if(length(startEnd)!=4 || var(sapply(startEnd, length))!=0 )
      stop("startEnd must have 3 vectors 'session', 'start', 'end' of equal length")
    sel = as.data.frame(startEnd,stringsAsFactors = F)
  } else stop("startEnd must be a list or dataframe")
  
  
  # if(!all.equal(names(sel),c('start', 'end')) ) stop("startEnd names must be 'start','end'")
  ref = 0
  sel[grepl("start|inizio",sel[,4],ignore.case = T),4] = start(x)[1]+start(x)[2]/frequency(x)
  sel[grepl("end|fine",sel[,4],ignore.case = T),4] = end(x)[1]+end(x)[2]/frequency(x)
  
  sel$start = timeMaster(sel$start, out="s")
  sel$end   = timeMaster(sel$end,   out="s")

  #check for artefact start lower than signal start
  if(any(sel$start < tss(x))){
    warning("artefacts times beginning before signal's start time were trimmed")
    sel[sel$start < tss(x),] = round(xstart(x))
  }
  #check for artefact end greater than signal end
  if(any(sel$end > tse(x))){
    warning("artefacts times ending after signal's end time were trimmed")
    sel[sel$end > tse(x),] = round(xend(x))
  }
  

  x$artefacts = sel[,c("start","end")]

  attributes(x)["filter"] = paste0(attr(x,"filter"), "--> artifacts set")
  x
}



############################################################################################################
############################################################################################################
############################################################################################################

#' @export
 signalDecimate = function (signal, newSampRate) {stop("this function has been deprecated. Use resample and signalFilter instead")}


#' resample
#' A simple wrapper for approx, used to decimate or upsample (by linear interpolation) a time series
#'
#' @param x A time-series object
#' @param newSampRate the new frequency
#' @param ... further options to be passed to approx
#'
#' @return a resampled time-series with the same start of the original time series.
#' @export
#'
#' @examples

resample = function (x, newSampRate, ...) {
  if(newSampRate == frequency(x)) stop("newSampRate and original signal frequency are identical.")
  if(newSampRate == 12) warning("by default, ts() assumes frequency Values of 4 and 12 to imply a quarterly and monthly series respectively (e.g.) in print methods.")
  q = frequency(x) / newSampRate  #ratio of old vs new sr
  ts(approx(seq_along(x),x, xout= seq(1,length(x),by=q), ... )$y, start=tss(x), frequency=newSampRate)
}
  





############################################################################################################
############################################################################################################
############################################################################################################
#### Moving average filters

#movAv has improved managing of burn-in and burn-out sequences
#if in doubt use this.

#' Moving average filter
#'
#' @param x A time-series or numeric vector
#' @param winSec Size of the moving window, in seconds
#' @param incSec Size of each window increment, in seconds. If NA a new window is
#' calculated for every sample.
#' @param remove If true, the function subtracts the rolling average from the
#' the original signal, acting as a high-pass filter. 
#' @param sampRate The frequency of x, i.e. samples per second
#'
#' @return A time-series or numeric vector of the same length of x.
#' @details The burn-in sequence has length of winSec/2 while the burn-out sequence might
#' be longer as not all combinations of winSec and incSec may perfectly cover any arbitrary x size. 
#' The burn-in and burn-out sequences are obtained through linear interpolation from their average value
#' to the first or last values of the moving average series.
#' The exact number of windows is given by ceiling((length(x)-win+1)/inc) where win and inc are
#' respectively the window and increment sizes in samples.
#' @export

movAv <- function(x, winSec, incSec = NA, remove=FALSE, sampRate=frequency(x) ) {

  win = winSec*sampRate
  win2 = round(win/2)
  inc = if(is.na(incSec)) 1 else incSec*sampRate
  len = length(x)
  n_win = ceiling((len-win+1)/inc)
  if(n_win<1) stop("During movAv routing the chosen window and increment led to zero windows.")
  winssamp = seq(1,by = inc, length.out = n_win)
  res = numeric(n_win)
  a = seq(1,by = inc, length.out = n_win)
  b = seq(0,by = inc, length.out = n_win) + win
  for(i in 1:n_win) {
    res[i] = sum(x[a[i]:b[i]], na.rm=T) /win
  }
  # if inc == 1 any window size will cover exactly all samples
  # because the last window will be (len-win+1):len
  # also, the resulting series has naturally the same sampling rate of the 
  # original series.
  
  # instead if inc is different than 1 the resulting series sampling rate
  # is equal to inc. The length of the resulting series is roughly len/inc
  # and needs to be interpolated back to the original sampling rate 
  
  if(inc>1){
    max_x = (n_win-1)*inc +1                    #starting value of the last window
    res = approx(x    = seq(1, max_x, by=inc), #starting sample of each window
                 y    = res,                    #average value of each window
                 xout = seq(1,max_x)            #x values of new series
    )$y
    
  }
  
  # Independently from inc, he resulting series will be (win-1) samples shorter than the original,
  # so it should be padded by win2 at the start and win2-1 end (plus eventual other missing samples
  # in the case of inc not multiple of len)
  # In this implementation the padding is a straight line from the mean of the starting and
  # ending parts of the original signal to the first/last calculated moving average point
  
  #pad initial half window
  startVal = mean(x[1:win2],na.rm=T)
  if(is.na(startVal)) startVal = res[1]
  res = c(seq(startVal,res[1],length.out = win2 ), res)
  
  #how many samples could not be estimated?
  miss = len - length(res) 
  #fill them with good values
  endVal = mean(x[(len-win2):len],na.rm=T)
  if(is.na(endVal)) endVal = res[length(res)]
  res = c(res,  seq(res[length(res)], endVal, length.out = miss))
  
  #all this done, x and res should ALWAYS have the same length
  if(length(res) != length(x)) warning("the original series and the moving average are of different lenghts, which is a bug")
  
  if(remove){
    res = x-res
  }
  
  if(is.ts(x)) res = ts(res,start=start(x),frequency = sampRate)
  else res
}


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
# stepCenter = function(a, winSec=60){
#   if(is.DyadStream(a)) x=a
#   
#   freq = frequency(a)
#   winSam = winSec*freq
#   n_win = ceiling((length(a)-winSam+1)/winSam)
#   resid = length(a)-n_win*winSam
#   for(i in 0:(n_win-1)){
#     al = (i*winSam +1)
#     bl = (i*winSam +winSam)
#     a[al:bl] = a[al:bl] - mean(a[al:bl])
#   }
#   if (length(resid)>0)
#     a[(length(a)-resid+1):length(a)] = a[(length(a)-resid+1):length(a)] - mean(a[(length(a)-resid+1):length(a)])
#   ifelse(is.DyadStream(x), cloneDyadStream(a,x) , a )
# }



#' FIR filter for physiological signals
#' This is a wrapper for signal::fir1, calculating decent values of filter order
#' based on engineering rules of thumb
#'
#' @param x a ts or DyadStream object
#' @param cut filter cut in Hz
#' @param type lowpass or highpass
#' @param NAsub signal::fir1 does not allow to have NAs. So they have to be substituted
#' @param attenDb attenuation in Db
#' @param burnSec head and tail of the filtered signals are bad. Burnsec specifies an amount of seconds
#'                over which the filtered signal is crossfaded with the original one.
#' @param plot logical. if TRUE plots frequency response, pass and stop bands
#' @param maxN maximum filter order
#'
#' @details The filter order is calculated automatically following "fred harris' rule of thumb"
#' and a maximum transition bandwidth of 1.2 times the cut frequency for lowpass (e.g. a cut at 10Hz will
#' achieve maximum attenuation at )
#' @export
#' @examples
#' Fs = 100; t = 5; samples = seq(0,t,len=Fs*t)
#' x = ts(sin(2*pi*3*samples + seq(-0.5,0.5,length.out = Fs*t), frequency = Fs)
#' x[1:50] = NA
#' xn = x+runif(1000,-0.5,0.5);
#' x1 = FIR(xn, "low", cut = 3,attenDb = 50,burnSec = 0, plot=TRUE)
#' plot(xn,col="grey");lines(x,col=3,lwd=2);lines(x1,col=2,lwd=3)

FIR = function(x, cut, type=c("low","high"), NAsub=NA, attenDb=50, burnSec = 0, maxN = 500, plot=F){
  # https://dsp.stackexchange.com/questions/37646/filter-order-rule-of-thumb 
  # https://www.allaboutcircuits.com/technical-articles/design-of-fir-filters-design-octave-matlab/
  # x = sin(1:1000/20)+seq(-0.5,0.5,length.out = 1000)
  # Fs = 1000;t = 2
  # samples = seq(0,t,len=Fs*t)
  # x = (sin(2*pi*100*samples) + sin(2*pi*120*samples)+ sin(2*pi*180*samples))/3
  # x2 = (sin(2*pi*180*samples))
  # x = ts(x, frequency = Fs);x2 = ts(x2, frequency = Fs)
  # plot(x,col="grey",xlim=c(1,1.1))
  # # x = x+runif(Fs*t,-0.5,0.5);
  # cut = 180

  x[is.na(x)]=NAsub
  type = match.arg(type, c("low","high"))
  
  max_band = cut*1.1
  delta_f = abs(max_band - cut) #abs is he same for high and low
  N = min(maxN,ceiling(attenDb * frequency(x) /(22*delta_f)))
  if(type=="high"){if(N%%2 != 0) N = N+1}
  if(type=="low" ){if(N%%2 == 0) N = N+1}
  wc = cut/(frequency(x)/2) #normalizza per nyquist freq.
  if(wc > 1) stop("Filter cut must be lower than Niquist frequency:",(frequency(x)/2))

  bf = signal::fir1(N,wc,type)
  if(plot){
    k = freqz(bf,Fc=frequency(x))
    freqz_plot(k)
  }
  xf = signal::filtfilt(filt = bf,x)

  #head and tail are bad, cross-fade the first and last n seconds
  if(burnSec>0){
    burn = burnSec*frequency(x)
    firstSamp = 1:(burn)
    lastSamp = (length(xf)-burn+1):length(xf) 
    
    burn1 =  x[firstSamp] * seq(1,0,length.out = burn) 
    burn2 =  x[lastSamp]  * seq(0,1,length.out = burn)
    
    #fadeout the filtered signal
    xf[firstSamp] = xf[firstSamp]*seq(0,1,length.out = burn) + burn1
    xf[lastSamp]  = xf[lastSamp] *seq(1,0,length.out = burn) + burn2
  
}
  #fir1 also changes the average value, recalibrate
  xf = xf-(mean(xf,na.rm=T) - mean(x,na.rm=T))
  if(length(xf)!=length(x)) warning("lenght")
  
  if(is.DyadStream(x))
    cloneAttr(x,ts(xf,frequency=frequency(x),start=start(x),end=end(x)))
  else 
    ts(xf,frequency=frequency(x),start=start(x))
}

x = ts(sin(1:1000/20)+seq(-0.5,0.5,length.out = 1000), frequency = 100)
# x[1:50] = NA
xn = x+runif(1000,-0.5,0.5);
x1 = FIR(xn, "low", cut = 3,attenDb = 50,burnSec = 0)
plot(xn,col="grey");lines(x,col=3,lwd=2);lines(x1,col=2,lwd=3)


