# ##PPG 2 HRV - 2024 version
# rm_all()
# load("A:/OneDrive - Università degli Studi di Padova/__Ricerche/2023_BIOFEEDBACK_PILOT2/PilotBF2_amico1.RData")
# library(DyadSync)
# 
# x = X = lr[[3]]$PPG$s1
# pk = peakFinder(x, mode = "peaks", FIR_args = list(type = "low", cut=3.0,NAsub = 0),
#                 correctionRangeSeconds = 1/1.5, minPeakAmplitude=0.01)
# 
# r = 0.65; f = 0.55; buffer  = 8;
# FIR_args = list(type = "low", cut=3.0,NAsub = 0);
# correctionRangeSeconds = 1/1.5; minPeakAmplitude=0.01;
# plotName = "test.svg"

#' @export
rejectFastPPG = function(pk,r=0.65, buffer = 8, plot=FALSE){
  pk2 = pk
  check = TRUE
  allRej = c()
  while(check){
    delta = diff(pk2$x)
    delta = c(NA,delta)
    toRej = c()
    for(i in seq_along(delta)){
      if(i>buffer){
        med = median(delta[(i-buffer):(i-1)],na.rm=T)
        if(delta[i]< med*r) {
          toRej = c(toRej, i)
          delta[i+1] = delta[i+1]+delta[i]
          delta[i] = NA
        }
      }
    }
    if(length(toRej) != 0){
      allRej = c(allRej, toRej)
      pk2$bool[pk2$samples[toRej]] = FALSE
      pk2$x = pk2$x[-toRej]
      pk2$y = pk2$y[-toRej]
      pk2$samples = pk2$samples[-toRej]
      pk2$time = pk2$time[-toRej]
      pk2$class = pk2$class[-toRej]
      pk2$amp = pk2$amp[-toRej]
      pk2$index = pk2$index[-toRej]
    } else {
      check = FALSE
    }
  }
  pk2$reject_samples = allRej

  if(!isFALSE(plot)){
    toPlot = allRej[c(3,diff(allRej))>2]
    for(j in sort(sample(seq_along(allRej),min(50,length(allRej))))){
      i = toPlot[j]
      plot(window(x, start =pk2$x[i-buffer*1.2],duration =buffer*1.3), main=paste("reject:",i))
      abline(v= pk$x,lty=3,col="#999999")
      abline(v= pk$x[c(i+1, i-1)],lty=1)
      abline(v= pk$x[allRej],lty=2, lwd=3, col=2)
    }
  }
  
  
  return(pk2)
}

#' @export
fillMissingPPG = function(x, pk, f = 0.75, buffer = 8, plot = FALSE){
  delta = diff(pk$x)
  delta = c(NA,delta)

  new_samples = c()
  new_i = c()
  for(i in seq_along(delta)){
    if(i>buffer){
      med = median(delta[(i-buffer):(i-1)],na.rm=T)
      if(delta[i] >= (med * f)*2){
        #how many points?
        n = floor(delta[i] / (med * f))
        for(j in 1:(n-1)){
          midpoint_s = pk$samples[i-1] + round((pk$samples[i]-pk$samples[i-1])/n)
          new_samples = c(new_samples, midpoint_s)
          new_i = c(new_i, i)
        }
        
      }
    }
  }
  
  pk2 = pk
  pk2$bool[new_samples] = TRUE
  pk2$samples = which(pk2$bool)
  pk2$x = c(pk2$x, x$x[new_samples])
  pk2$y = c(pk2$y, x$y[new_samples])
  pk2$y = pk2$y[order(pk2$x)]
  pk2$class = c(pk2$class, rep("p_filled",times = length(new_samples)))
  pk2$class = pk2$class[order(pk2$x)]
  pk2$amp = c(pk2$amp, NA)
  pk2$amp = pk2$amp[order(pk2$x)]
  pk2$x = sort(pk2$x)
  pk2$time =  pk2$x
  pk2$index = c(pk2$index, pk2$index[length(pk2$index)]+cumsum(rep(2,length(new_samples))))
  
  if(!isFALSE(plot)){
    for(j in seq_along(new_samples)){
      i = new_i[j]
      plot(window(x, start =pk2$x[i-buffer*1.2],duration =buffer*1.3), main=paste("fill:",i))
      abline(v= pk$x,lty=3,col="#999999")
      abline(v= pk2$x[c(i+1, i-1)],lty=1)
      abline(v= pk2$x[i], col=3,lty=2, lwd=3)
      points(x$x[pk2$bool], x$y[pk2$bool], col=4)
      points(x$x[pk2$samples[i]], x$y[pk2$samples[i]], col=2,cex=2,pch=15)
      points(pk2$x[i], pk2$y[i], col=3,cex=2,pch=15)
      
    }
  }
  return(pk2)
}


#' Extract beats from PPG signal
#' @description
#' This functions employs a low pass filter to attenuate non-systolic peaks,
#' then rejects peaks that are too fast compared to the median rate of a buffer.
#' Finally it interpolates too large gaps.
#' 
#' @param x a rats time-series representing the original PPG signal
#' @param r the rejection factor, expressed as a weight of the buffer's median R-R interval 
#' @param f the fill factor, expressed as a weight of the buffer's median R-R interval 
#' @param buffer integer. the rolling buffer, expressed in number of beats
#' @param FIR_args a named list of arguments to be passed to \link[DyadSignal]{FIR}) filtering function.
#' 'cut' and 'type' must be present, representing . If empty, the filtering is disabled.
#' @param correctionRangeSeconds the half range in which the real maximum/minimum
#' value should be searched, in seconds.
#' around the derivative shift. Should be less then the periodicity of the
#' signal.  0.5 is good for skin conductance.
#' @param minPeakAmplitude the minimum delta from valley to peak for detection.
#' Skin conductance traditionally uses 0.05uS, or 0.03uS
#' @param plotName an optional filename to save a svg representation of the beats extraction
#'
#' @returns
#' a list of:
#' \describe{
#'   \item{bool}{a logical vector of the same length of x with TRUE
#'    value corresponding to a match
#'    }
#'   \item{samples}{the position of beats relatively to x's index}
#'   \item{x}{temporal coordinates of the beats (if x has a frequency attribute)}
#'   \item{y}{The value of x correspoding to the detected beats}
#'   \item{type}{a character vector defining for each beats if its a detected peak
#'               ('p') or an interpolated one ('p_filled')}
#'   \item{amp}{the through-to-peak amplitudes of the detected beats}
#'   \item{index}{The absolute position of the beats along the signal. The values
#'      are independent from the 'mode' argument, allowing interaction between
#'      different calls of the function.}
#'   \item{reject_samples}{the position of rejected beats relatively to x's index}
#'   
#'   
#' }
#' @export
#'
#' @examples
PPG2beats = function(x, r = 0.65, f = 0.75, buffer  = 8,
                   FIR_args = list(type = "low", cut=3.0,NAsub = 0),
                   correctionRangeSeconds = 1/1.5, minPeakAmplitude=0.01,
                   plotName = NULL, ...
                   ){
  UseMethod("PPG2beats", x)
}

#' @export
PPG2beats.DyadExperiment = function(x, r = 0.65, f = 0.75, buffer  = 8,
                          FIR_args = list(type = "low", cut=3.0,NAsub = 0),
                          correctionRangeSeconds = 1/1.5, minPeakAmplitude=0.01,
                          plotName = NULL, ...
                          ) {
  
  for(nSession in seq_along(x)) {
    ppg_i = grep("ppg", names(x[[nSession]]), ignore.case = TRUE)
    prog(nSession,length(x))
    
    # print(nSession)
    
    if(length(ppg_i) != 1) stop("None or multiple DyadSignals contained PPG in the name")
    
    if(!is.null(plotName)) dir.create(plotName, showWarnings = FALSE, recursive = TRUE)
    plotNameS1  = if(!is.null(plotName))  file.path(plotName, paste0(UID(x[[nSession]]),"_","s1" )) else NULL
    plotNameS2  = if(!is.null(plotName))  file.path(plotName, paste0(UID(x[[nSession]]),"_","s2" )) else NULL
    
    tryCatch({

      
      b1 = PPG2beats(x[[nSession]][[ppg_i]]$s1, r = r, f = f, buffer  = buffer,
                     FIR_args = FIR_args,
                     correctionRangeSeconds = correctionRangeSeconds, minPeakAmplitude=minPeakAmplitude,
                     plotName = plotNameS1, ...
      )
    }, error = function(e) {
      b1 = NA
      warning("The following issue was found for s1 in session:", UID(x[[nSession]]),"\r\n",e,"\r\n" )
    })
    
    tryCatch({
      b2 = PPG2beats(x[[nSession]][[ppg_i]]$s2, r = r, f = f, buffer  = buffer,
                     FIR_args = FIR_args,
                     correctionRangeSeconds = correctionRangeSeconds, minPeakAmplitude=minPeakAmplitude,
                     plotName = plotNameS2, ...
      )
    }, error = function(e) {
      b2 = NA
      warning("The following issue was found for s2 in session:", UID(x[[nSession]]),"\r\n",e,"\r\n" )
    })
    
    
    x[[nSession]][[ppg_i]][["beats"]] = list(s1 = b1, s2 = b2)
  }
  return(x)
}

#' @export
PPG2beats.rats = function(x, r = 0.65, f = 0.75, buffer  = 8,
                     FIR_args = list(type = "low", cut = 3.0,NAsub = 0),
                     correctionRangeSeconds = 1/1.5, minPeakAmplitude=0.01,
                     plotName = NULL
                     ) {
  buffer = trunc(buffer)
  pk = peakFinder(x, mode = "peaks", FIR_args = FIR_args,
                  correctionRangeSeconds = correctionRangeSeconds, minPeakAmplitude = minPeakAmplitude)
  step1 = rejectFastPPG( pk, r = r, buffer = buffer, plot = FALSE)
  step2 = fillMissingPPG(x, step1, f = f, buffer = buffer, plot = FALSE)

  if (!is.null(plotName)) {
    svg(filename = plotName, height = 10, width = duration(x)/2)
    plot(x,xaxs = "i")
    abline(v = step2$x,lwd = 2)
    abline(v = pk$x[step1$reject_samples], col = 2, lwd = 2, lty = 2)
    abline(v = step2$x[step2$class == "p_filled"], col = 3, lwd = 2, lty = 2)
    dev.off()
  }

  return(step2)
}

beats2HR = function(experiment , signalName = "HR", ...){
  for (nSession in seq_along(experiment )) {
    ppg_i = grep("ppg", names(experiment [[nSession]]), ignore.case = TRUE)
    prog(nSession,length(experiment ))
    if (length(ppg_i) != 1) stop("None or multiple DyadSignals contained PPG in the name")
    b1 = experiment [[nSession]][[ppg_i]][["beats"]]$s1
    b2 = experiment [[nSession]][[ppg_i]][["beats"]]$s2
    
    # @HACK there should not be repeating experiment  values in peaks
    # yet ....
    #rle removes duplicate values
    n = 10
    # lim = 200
    x1 = rle(b1$x)$values  
    rr1 = diff(x1)
    rr1f = filter(rr1, rep(1 / n, n), sides = 2)
    # plot(60/rr1,t="l",xlim=c(1, lim));lines(60/rr1f, col=2)
    x2 = rle(b2$x)$values  
    rr2 = diff(x2)
    rr2f = filter(rr2, rep(1 / n, n), sides = 2)
    # plot(60/rr2,t="l",xlim=c(1, lim));lines(60/rr2f, col=2)
    
    plot(x2[2:length(x2)],rr2f)
    lines(RR2,col=2)
    RR1 = rats(approx(xout=experiment[[nSession]][[ppg_i]]$s1$x, x = x1[2:length(x1)], y = rr1f )$y,
               start = start(experiment[[nSession]][[ppg_i]]$s1),
               frequency = frequency(experiment[[nSession]][[ppg_i]]$s1),
               timeUnit = timeUnit(experiment[[nSession]][[ppg_i]]$s1),
               unit="beats/min")
    RR2 = rats(approx(xout=experiment[[nSession]][[ppg_i]]$s2$x, x = x2[2:length(x2)], y = rr2f )$y,
               start = start(experiment[[nSession]][[ppg_i]]$s2),
               frequency = frequency(experiment[[nSession]][[ppg_i]]$s2),
               timeUnit = timeUnit(experiment[[nSession]][[ppg_i]]$s2),
               unit="beats/min")
    newSignal = DyadSignal(name= signalName,s1 = 60/RR1, s2 = 60/RR2, SR = frequency(experiment[[nSession]][[ppg_i]]$s1),
                s1Name = s1Name(experiment[[nSession]]),
                s2Name = s2Name(experiment[[nSession]]),
                sessionId = sessionId(experiment[[nSession]]), dyadId = dyadId(experiment[[nSession]]), groupId=groupId(experiment[[nSession]]))
    experiment[[nSession]][[signalName]] = newSignal
    ###
    # HRV computation
    #https://www.kubios.com/blog/hrv-analysis-methods/
    HRV_time = function(x){
      #NB x must be a RR sequence
      # x = RR2
      # x = rr2f
      #to milliseconds

      x = x * 1000
      time_domain = list(
        RMSSD = sqrt(sum(diff(x)^2,na.rm=T)/length(x-1)),
        meanHR = 60/mean(x,na.rm=T),
        SDNN =  sqrt(sum((x-mean(x, na.rm=T))^2, na.rm=T )/length(x-1))
        )
      time_domain
    }

    
    
    HRV_time(rr2f)
    hd = RHRV::CreateHRVData()
    hd = RHRV::SetVerbose(hd,TRUE)
    
    

    dat = lr2[[2]]$PPG$s2
    test = PPG2beats(dat,plotName = "testPPG.svg")
    plot(window(dat,start= 600, duration = 30))
    hd = RHRV::LoadBeatVector(hd, test$x)
    hd = RHRV::BuildNIHR(hd)
    hdf = RHRV::FilterNIHR(hd)
    rr2f = filter(diff(test$x), rep(1 / n, n), sides = 2)
    
    abline(v = test$x ,lwd=2)
    abline(v = hdf$Beat$Time, col=2,lty=2,lwd=2)
    abline(v = cumsum(na.omit(rr2f )), col=4,lty=2)
    
    
    hd  = RHRV::CreateTimeAnalysis(hd,size=300,interval = 7.8125)
    HRV_time(diff(test$x))
    hdf = RHRV::CreateTimeAnalysis(hdf,size=300,interval = 7.8125)
    
    
    RRDiffs = diff(hd$Beat$RR)
    sqrt(mean(RRDiffs^2))
    RRDiffs = diff(hdf$Beat$RR)
    sqrt(mean(RRDiffs^2))
    RRDiffs = diff(rr2*1000)
    sqrt(mean(RRDiffs^2))
    RRDiffs = diff(na.omit(rr2f)*1000)
    sqrt(mean(RRDiffs^2))
    
    HRV_time(hd$Beat$niHR)
    HRV_time(hdf$Beat$niHR)
    
    
    hd = RHRV::InterpolateNIHR (hd, freqhr = 10)
    
    hd = RHRV::CreateFreqAnalysis(hd)
    hd$TimeAnalysis
    
    
    
    lapply(HRV_time(rr2f), \(x){x*1000})
    HRV_freq = function(x) {
      x2 = movAv(x, winSec = 120, incSec=1/frequency(x),remove = T, SR = frequency)
      
    }
      #detrend 
      plot(x2)
      }
    
  
}

#' 

#' 
#'     
#'     hr1 = rats(approx(xout=x$s1$x, x = b1$x[2:length(b1$x)], y = diff(b1$x) )$y,
#'                start = start(x$s1), frequency = frequency(x$s1), timeUnit = timeUnit(x$s1),
#'                unit="R-R interval")
#'     hr2 = rats(approx(xout=x$s1$x, x = b1$x[2:length(b1$x)], y = diff(b1$x) )$y,
#'                start = start(x$s1), frequency = frequency(x$s1), timeUnit = timeUnit(x$s1),
#'                unit="R-R interval")
#'     
#'     newSignal = DyadSignal(name= "heart",s1 = s1raw, s2 = s2raw, SR = jSampRate,
#'                            s1Name= paste0(comboList[[j]][i,"ranx_session"], if(comboList[[j]][i,"ranx_role"] =="s1") olds1name else olds2name ),
#'                            s2Name =paste0(comboList[[j]][i,"rany_session"], if(comboList[[j]][i,"rany_role"] =="s1") olds1name else olds2name ),
#'                            sessionId = i, dyadId = "random", groupId="random")
#'     
#'     signalList[[j]] = newSignal
#'     
#'     
#'     prog(nSession, length(x))
#'     
#'     
#'     
#'     x[[nSession]] = epochSeries(x[[nSession]], signal=signal, sync=sync, series=series,
#'                                 diary=diary, category=category, mergeEpochs=mergeEpochs,
#'                                 artefact.rm=artefact.rm, shift_start=shift_start, shift_end=shift_end, target=target,
#'                                 summaryFUN = summaryFUN, ...)
#'   }
#'   
#'   return(x)
#'   
#'   
#' }
#' 
