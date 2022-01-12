##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################                                   
## DyadUtils.R
# good functions for good work
#
############################################################################################################
## Changelog
# v1.2 - small (but significant) edits to timeMaster()
# v1.1 - partially compatible wiht dyadClass v1.3 (no more sessions)
# v1.0 - fully 'stream 1.1' compliant
#
############################################################################################################
## ToDo
#
#
############################################################################################################ 
## Credits
# Author: Johann R. Kleinbub
# Contact: johann.kleinbub@gmail.com
############################################################################################################ 



############################################################################################################
############################################################################################################
############################################################################################################
## FUNCTIONS FOR EXPERIMENTS


##
#' THIS IS OLD
#' USE signalFilter instead
#' applies a function on each signal of an experiment.
#' @param experiment 
#' @param signals either the string "all" or a vector of signal names.
#' If 'all' the function will be applied on all DyadSignal objects. To apply on epochs,
#' and other objects, the names must be specified
#' @param FUN a function passed to a lapply call on a list of DyadSignals
#' @param ... further arguments passed to lapply
#'
#' @details the FUN can be in the generic form function(x, ...){ } where x is each DyadSignal object
#' among the selected signals, in each session.
#' @export
#'
#' @examples
expApply = function(experiment, signals="all", FUN, ...){
  warning("expApply just applies a function on each session. Try signalFilter!")
  #signals can either be "all" or a vector of signal names. es: c("PPG","SC").
  if(!is(experiment,"DyadExperiment")) stop("Only objects of class DyadExperiment can be processed by this function")
  fun = match.fun(FUN)
  experiment2 = Map(function (session,nSession){
    ##debug
    #     session = lr$sessions[[1]]
    #     signals= c("SC","ASD","PPG")
    
    if(length(signals)==1 && signals=="all") {
      sigs = names(session)[sapply(session,is.DyadSignal)]
      } else sigs = signals
    session[names(session) %in% sigs] = lapply(session[names(session) %in% sigs],  fun, ...)
    prog(nSession,length(experiment))
    return(session)
  },experiment,seq_along(experiment))
  attributes(experiment2)=attributes(experiment)
  return(experiment2)
}

# DELETEME
#' @export
experimentMerge = signalMerge = sessionMerge = function(...){
  stop("These functions are obsolete. Please use c(...)")
}


##
#' Title
#' This function allows to subset only specific signals from an experiment
#' @param x a DyadExperiment or DyadSession object 
#' @param signals a vector of signal names. Es: c("SC", "HRV", "PACS")
#'
#' @return
#' @export
#'
#' @examples
selectSignals = function(x, signals){
  UseMethod("selectSignals",x)
}
#' @export
selectSignals.DyadSession = function(session,signals) {
  session[!names(session) %in% signals] = NULL
  session
}
#' @export
selectSignals.DyadExperiment = function(experiment,signals) {
  res = lapply(experiment, selectSignals, signals)
  DyadExperiment(name(experiment),res)
}




# ccfQuantile = function (EXP, signal="SC",lag="bestCCF",bySession=T,sessionFUN = "mean"){
#   FUN = match.fun(sessionFUN)
#   resList = lapply(EXP, function(ses){ses[[signal]]$ccf$ccfmat[[lag]]})
#   res = do.call("rbind",resList)
#   if(bySession) res =apply(res,1,FUN)
#   quantile(res)
# }




#' Title
#'
#' @param a 
#' @param b 
#'
#' @return
#' @export

cohend = function(a,b,na.rm = TRUE){
  if(na.rm==T){
    a = a[!is.na(a)]
    b = b[!is.na(b)]
  }
  pool = sqrt( (  (length(a)-1)*sd(a)^2 +  (length(b)-1)*sd(b)^2 )/(length(a)+ length(b)-2 ) )
  (mean(a)-mean(b))/pool
}

#' pvalFormat
#'
#' @param x 
#'
#' @return
#' @export

pvalFormat= function(x){
  if(x>= 0.0001) {
    if(x<0.001){
      x = "< 0.001 ***"
    } else if(x<0.01){
      x = paste(round(x,3),"**")
    } else if(x<0.05){
      x = paste(round(x,3),"*")
    } 
  } else x = "< 0.0001 ****"
  x
}





############################################################################################################
############################################################################################################
############################################################################################################
## GLOBAL TOOLS

#' Extends tolower toupper with title case translation
#' @export
totitle <- function(x) {
  sapply(x, function(k){
    s <- strsplit(k, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
          sep = "", collapse = " ")
  },USE.NAMES =F)

}

#' @export
cat0 = function(...) {cat(..., sep="")}

#' @export
lead0 = function(x, width = 2){
  formatC(x,width = width, format = "d", flag = "0")
}

#' rangeRescale
#' this extremely useful function rescales a numeric vector to a new range through linear interpolation
#' @param x 
#' @param newA,newB  Any two points (typically min and max) of the new scale
#' @param oldA,oldB  The corresponding two points of the original scale.
#' If oldA or oldB are missing the min or max of x are used instead
#' @param pres.sign logical. if TRUE the signs of the original data are preserved. newA and newB must be opposites (e.g. -1,1)
#' @details 
#' @export
#'
#' @examples
rangeRescale <- function(x, newA, newB, oldA =min(x, na.rm=T), oldB = max(x, na.rm=T), pres.signs=FALSE){
  #newA e newB indicano il minimo e il massimo della nuova scala
  #
  #se oldA e oldB mancano, vengono usati il minimo e il massimo del campione
  # if(any(x>oldB, na.rm=T) || any(x<oldA, na.rm=T)) stop ("Value found outside oldB and ymin boundaries. oldB and oldA should be equal or larger than the data range.")
  if(pres.signs){
    # if((!missing(oldA) && !missing(oldB)) || (oldA!=-max(abs(x)) || oldB != max(abs(x)) ) )
    #   stop("Either x")
    if( newA != -newB )
      stop("with pres.sign = TRUE, newA and newB should be opposites (e.g. -1 and 1")
    mightyMax = max(abs(oldB),abs(oldA)) #così centra qualsiasi range?
    oldA = -mightyMax
    oldB = mightyMax
  }
  (newB-newA) * (
    (x - oldA)  / (oldB - oldA )
  ) + newA
}


prog = function(i,n,step=50){
  progStep = c(1,round((n/step)*2:(step-1)),n)
  progStep[progStep==0] = 1
  s = sum(i == progStep)
  if(s) {
    if(i==1) cat0(rep('.',50),"|100%\r\n")
    cat0(rep('.',s))}
  if(i==n) cat0("|Done ;)\r\n")
  
}

#This function checks if a package is present in R library
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

## 
#' cbind allowing to bind unequal columns to a same data.frame padding NAs to the end of the shorter
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
unequalCbind = function(...) {
  dots <- list(...)
  #debug
  #dots = list(my.orig,my.ccf)
  #
  dots = dots[!sapply(dots,is.null)]
  # print(str(dots,max.level=2))
  
  #print(str(dots))
  if(length(dots)>1){
    dots = lapply(dots, function(x){if(!is.data.frame(x)) data.frame(x) else x })
    dotsNames = unlist(sapply(dots,colnames))
    #print(dotsNames)
    maxlen = max(sapply(dots, nrow))
    fdots = lapply(dots, function(x){
      #deb
      #x= dots[[2]]
      #rm(x,y,fdots,dots,pad,maxlen)
      if(nrow(x)<maxlen){
        pad = maxlen - nrow(x) 
        y = data.frame(matrix(rep(NA,ncol(x)*pad),ncol = ncol(x))) 
        colnames(y) = colnames(x)
        rownames(y) = paste0("NA",1:pad)
        rbind(x,y)
      } else x
    })
    res = data.frame(do.call("cbind",fdots))
    colnames(res) = dotsNames
    #print(colnames(res))
    return(res)
  } else return(dots[[1]])
}


#' timeMaster
#' timeMaster allows to:
#' transform time from different formats to different formats.
#' add amounts of time to baseTime expressed in different formats.
#' This is verified. And useful. For reasons. :-D
#'
#' @param baseTime 
#' @param out 
#' @param add 
#' @param baseSep 
#'
#' @return
#' @export
#'
#' @examples
timeMaster = function(baseTime, out=c("auto", "hour", "min","sec"), add=0, baseSep = "[\\.,:,;,\\,',\",\\-]"){
  #baseTime and add can either be  integers of seconds or a time string in the format h:m:s, m:s, or s, with or without leading zeroes
  #output forces the sum to be reported either as string h:m:s or m:s or as a integer of seconds. auto keeps the 'baseTime' format.
  out = match.arg(out)
  if(is.numeric(baseTime) && sum(baseTime%%1)>0 ) {
    baseTime=as.character(baseTime)
    warning("The . was considered as a ':' in the format mm:ss. timeMaster does not support fractional times yet")
  }
  if(length(baseTime)>1)
  {
    ## NB questa è la linea classica, funzionantissima, tranne nel caso di un data.frame di una riga.
    SIMPLIFY = if(is.data.frame(baseTime)) FALSE else TRUE
    sapply(baseTime,timeMaster,out,add,baseSep,USE.NAMES = F)
    
    ## Questo è il nuovo approccio. Potrebbe sfasciare tutto
    # res = baseTime
    # for(i in seq_along(baseTime)){
    #   res[[i]] <- timeMaster(baseTime[[i]],out,add,baseSep)  
    #   print
    # }
    # res
    
    ## oppure
    # if(is.list(baseTime)) res = vector("list", length(baseTime)) 
    # else res = numeric(length(baseTime))
    # for(i in seq_along(baseTime)){
    #   res[[i]] <- timeMaster(baseTime[[i]],out,add,baseSep) 
    #   # cat("\r\n",str(res))
    # }
    # names(res) = names(baseTime)
    # # class(res) = class(baseTime)
    # res
    # 
    ## elimina fino a qui in caso.
    
  } else {
    #da qui baseTime è contenente un tempo singolo, non un vettore
    if(is.character(baseTime)){
      #è negativo?
      if(substr(baseTime,1,1)=="-"){
        negative = T
        baseTime = substring(baseTime,2)
      } else negative = F
      baseTime = strsplit(baseTime, split=baseSep)
      if      (length(baseTime[[1]])==1) auto = "sec"
      else if (length(baseTime[[1]])==2) auto = "min"
      else if (length(baseTime[[1]])==3) auto = "hour"
      else stop ("baseTime format not recognized. It should either be \"min:sec\" or \"hour:min:sec\" or an integer of seconds")
      sapply(unlist(baseTime), function(k){if(k=='') stop("baseTime contains an empty cell. Please check your separators")})
      #transform to seconds
      baseTime = unlist(lapply(baseTime, function(x){
        if(length(x)==1)
          as.numeric(x[1])
        else if(length(x)==2)
          as.numeric(x[1])*60 + as.numeric(x[2])
        else if(length(x==3))
          as.numeric(x[1])*3600 + as.numeric(x[2])*60 + as.numeric(x[3])
        else stop("Time format must either be \"min:sec\" or \"hour:min:sec\"")
      } ))
      if(negative) baseTime = -baseTime
    } else {
      auto = "sec"
    }
    if(is.character(add)) add = timeMaster(add,out="sec")
    x = baseTime + add #that's the final amount in seconds
    
    if(out=="auto") out=auto
    if(out == "sec") {
      return(x)
    } else if(out %in% c("min", "hour")){
      if (x<0) { x = abs(x); negative = T} else negative =F
      mins = floor(x/60)
      secs = (x-mins*60)
      hours = trunc(mins/60)
      mins_h = mins - hours*60
      if(out == "min")
        return(paste0(ifelse(negative,"-",""), lead0(mins),":", lead0(secs)))
      else if (out =="hour"){
        return(paste0(ifelse(negative,"-",""),lead0(hours),":",lead0(mins_h),":", lead0(secs)))
      }
    } else stop("timeMaster failed in a weird way!")
  }
}

#' @export
merge.list = function(x,y) {
  # l = list(...)
  xNames = names(x)
  yNames = names(y)
  over = names(x)[which(names(x) %in% names(y))]
  for(n in over){if(x[[n]]!=y[[n]]) stop("Can't merge list with different values for the same tag:\r\n",n,": ",x[[n]]," != ",y[[n]])}
}

#' Calculate number of windows in a time-series with a given duration
#'
#' @param x either the duration in seconds of a time-series, or a ts object
#' @param winSec the width of each window, in seconds
#' @param incSec the amount of increment between each window (incSec == winSec gives non-overlapping windows)
#' @param sampRate the number of samples per second (i.e. frequency)
#' @param return either "number", which returns the number of windows, or "all" which returns a data.frame with the start and end of each window.
#' @param verbose should all the windows be printed?
#' @param flex if true the first and last windows are stretched so to have exactly length(x)/inc resulting windows
#'
#' @return
#' @export
#'
#' @examples
nwin = function(x, winSec, incSec, sampRate=frequency(x), return=c("number", "all"), flex=FALSE, verbose = FALSE) {
  if(is.ts(x)){
    len = length(x)
  } else if(length(x)>1 || x%%1 != 0){
      stop ("in nwin() x must be either a ts object or a single integer number")
  } else {
      if(missing(sampRate)) warning("in nwin() sampRate was coerced to 1")
      len = x *sampRate
  }



  # sampRate = 10
  inc=incSec *sampRate;
  win=winSec*sampRate;


  #la vera santa formula:
  n_win = ceiling((len-win+1)/inc)

  return=match.arg(return)

  all_wins = 1:n_win
  start = (all_wins-1)*inc +1
  end = start + win -1
  mid = start + win/2

  # if(flex) {
  #   init =
  # }



  if(verbose){
    digits = nchar(trunc(abs(len)))

    cat("\r\n Number of windows:",n_win)
    cat("\r\n Number of samples used:",(n_win-1)*inc+win)
    cat("\r\n Number of seconds used:",((n_win-1)*inc+win)/sampRate)
    cat0("\r\n Proportion of len used: ",((n_win-1)*inc+win)/len*100,"%\r\n")

    actual = 0;iterations = 0;
    while ( actual+win<=len) {
      iterations = iterations+1
      cat0('\r\nW',lead0(iterations, nchar(n_win))," | Samples: ", lead0(actual+1,digits),' - ',lead0(actual+win,digits),
           " | Seconds: ", lead0(actual+1/sampRate,digits),' - ',lead0((actual+win)/sampRate,digits))
      actual = actual+inc;
      #cat(actual+win,"\n\r")
    }
  }

  if(return=="number")  return(n_win)
  else if (return =="all") return(data.frame(window=all_wins, start=start,end=end, mid=mid))
}


#' Interpolate by window
#' Interpolate windowed data to original samplerate. 
#' useful to overlay computed indexes on original timeseries
#'
#' @param windowsList
#' @param winSec 
#' @param incSec 
#' @param sampRate 
#'
#' @return
#' @export
#'
#' @examples
winInter = function(windowsList, winSec, incSec, sampRate){
  warning("This function is a mess, please ask some developer to refactor it!")
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

#' Approx by window
#' when you have a value for each window (even overlapping ones) and you
#' want to transform the values to a time series with the original sampling rate
#'
#' @param a a numerical vector, resulting from a windowing
#' @param winSec the size of the window used to calculate a, in seconds
#' @param incSec the increment of each window used to calculate a, in seconds
#' @param sampRate the number of samples per second of the destination time-series
#' @param midPoints if true, the value are considered to be at the center of the window, and the result is padded right by half window
#'
#' @return
#' @export
#'
#' @examples

# approxByWin = function(a, winSec, incSec, sampRate=frequency(x), midPoints=T){
#   # a = res
#   # winSec=4
#   # incSec = 1
#   # sampRate  = 100
#   warning("I think this function has issues!")
# 
#   inc=incSec*sampRate
#   win= winSec*sampRate
# 
#   halfWin = ceiling(win/2)
#   if(midPoints)
#     x_vals = seq(1,(length(a)*inc), by=inc)+halfWin #windowed datapoints must be at the center of the window
#   # else
#   #   x_vals = seq(1,(length(a)*inc), by=inc) #windowed datapoints must be at the center of the window
# 
#   out = approx(x    = x_vals,
#                y    = a, #the exact calculated values
#                xout = seq(1,(length(a)-1)*inc + win,by=1) #all samples covered by the windows
#   )$y
#   #replicate first/last good value in first half windows
#   out[1:halfWin-1] = out[halfWin]
#   out[(length(out)-halfWin+1):length(out)] = out[(length(out)-halfWin)]
#   return(out)
# }


