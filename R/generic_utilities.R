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

#' cat without spaces
#' @param ... 
#'
#' @export
cat0 = function(...) {cat(..., sep="")}

#' Formats number with leading zeros
#' @param x 
#'
#' @param width number of leading zeros
#' @param digits number of decimal digits
#'
#' @export
lead0 = function(x, width = 2, digits=0){
  formatC(x,width = width, format = "f", digits=digits, flag = "0")
}

#' rangeRescale
#' this function rescales a numeric vector to a new range through linear interpolation
#' @param x 
#' @param newMin,newMax  Any two points (typically min and max) of the new scale
#' @param oldMin,oldMax  The corresponding two points of the original scale.
#' If oldMin or oldMax are missing the min or max of x are used instead
#' @param pres.sign logical. if TRUE the signs of the original data are preserved. newMin and newMax must be opposites (e.g. -1,1)
#' @details 
#' @export
#'
#' @examples
rangeRescale <- function(x, newMin, newMax, oldMin =min(x, na.rm=T), oldMax = max(x, na.rm=T), pres.signs=FALSE){
  #newMin e newMax indicano il minimo e il massimo della nuova scala
  #
  #se oldMin e oldMax mancano, vengono usati il minimo e il massimo del campione
  # if(any(x>oldMax, na.rm=T) || any(x<oldMin, na.rm=T)) stop ("Value found outside oldMax and ymin boundaries. oldMax and oldMin should be equal or larger than the data range.")
  if(pres.signs){
    # if((!missing(oldMin) && !missing(oldMax)) || (oldMin!=-max(abs(x)) || oldMax != max(abs(x)) ) )
    #   stop("Either x")
    if( newMin != -newMax )
      stop("with pres.sign = TRUE, newMin and newMax should be opposites (e.g. -1 and 1")
    mightyMax = max(abs(oldMax),abs(oldMin)) #così centra qualsiasi range?
    oldMin = -mightyMax
    oldMax = mightyMax
  }
  return((newMax-newMin) * ((x - oldMin)  / (oldMax - oldMin )) + newMin)
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
#' @param x either the duration in seconds of a time-series, or a rats or ts object
#' @param winSec the width of each window, in seconds
#' @param incSec the amount of increment between each window (incSec == winSec gives non-overlapping windows)
#' @param SR the number of samples per second (i.e. frequency)
#' @param return either "number", which returns the number of windows, or "all" which returns a data.frame with the start and end of each window.
#' @param verbose should all the windows be printed?
#' @param flex if true the first and last windows are stretched so to have exactly length(x)/inc resulting windows
#'
#' @return
#' @export
#'
#' @examples
nwin = function(x, WIN, INC, flex=FALSE, SR=frequency(x), return=c("number", "all"),
                verbose = FALSE) {
  # DEBUG
  # stop("debug nwin")
  # x = lr10$all_CC_1$SC$s1
  # x = 1:100
  # SR = frequency(x)
  # WIN = 13
  # INC= 2
  # flex=T
  # return="all"
  # verbose=T
  # nwin(x=rats(1:1000, frequency = 0.16, start=17.33),WIN=21.17,INC=3.14,flex=T,return="all")
  #######
  
  if(length(x)>1){
    len = length(x)
  } else {
    if(x%%1 != 0) stop = "in nwin, x must be integer"
    len = x*SR
    x = 1:len
  }
  tstart = start(x)
  if(length(tstart) > 1){
    if(!is.null(attr(x, "tsp"))){ 
      tstart = tsp(x)[1]
      } else tstart = 1
  }

  # SR = 10
  inc=INC * SR;
  win=WIN*SR;
  if(inc<1 || win<1) stop("The data sampling rate is not sufficient to represent such small windows or increment")
  if(inc%%1!=0 || win%%1!=0) warning("Windows size or increment are not multiple of sampling rate. Calculation was approximated")
  inc = round(inc)
  win = round(win)
  if(win == 0 || inc==0) stop("Windows size and increment must be of at least 1 sample")
  win2 = win/2;
  
  #la vera santa formula:
  n_win = ceiling((len-win+1)/inc)

  return=match.arg(return)

  all_wins = 1:n_win
  start = (all_wins-1)*inc +1
  end = start + win -1
  mid = start + win2
  FLX = rep(FALSE, n_win)

  if(flex && win > inc){ 
    #flex is a system to have dynamically smaller windows to have exactly 1 value every inc
    #even if the windows size is great:
    
    #' \._./\.^.-^\._./\.^.-^\._./\.^.-^\._./\.^.-^
    #' |---------------o---------------|
    #'     |---------------o---------------|
    #'         |---------------o---------------|
    #'             |---------------o---------------|
    #'             
    #' Here you need to have at least half window before getting reliable values
    #' Instead with flex:
    #' \._./\.^.-^\._./\.^.-^\._./\.^.-^\._./\.^.-^
    #' o---------------|                  <---FLEX WINDOW 1
    #' |---o---------------|              <---FLEX WINDOW 2
    #' |-------o---------------|          <---FLEX WINDOW 3
    #' |-----------o---------------|      <---FLEX WINDOW 4
    #' 
    #' |---------------o---------------|  <---FIRST COMPLETE WINDOW
    #'     |---------------o---------------|
    #'         |---------------o---------------|
    #'             |---------------o---------------|
    #' 
    
    #initial flexes
    nif = floor((mid[1]) / inc) #number of flex windows
    if(nif > 0){
      fistart = rep(1,nif)
      fiend = seq(end[1]-inc*nif ,end[1]-1, by=inc)
      fimid = (fiend - win2)+1
      fiflex =rep(TRUE, nif)
      
      #@HACK non sono capace di trovare il numero giusto di nif
      #quindi correggo a mano.... 
      td = which(fimid<1)
      if(length(td)>0){
        fistart = fistart[-td]
        fimid   = fimid[-td]
        fiend   = fiend[-td]
        fiflex  = fiflex[-td]
      }

      

      #final flexes using also the unused data at the end
      end[n_win]
      start[n_win]
      
      ffstart = seq(start[n_win]+inc, length(x)-win2, by=inc)
      ffstart =cumsum(rep(inc,))
      
      nff = length(ffstart)
      ffend = rep(length(x), nff)
      ffmid = ffstart + win2
      ffflex = rep(TRUE, nff)
      
      start = c(fistart, start, ffstart)
      mid   = c(fimid,   mid  , ffmid)
      end   = c(fiend,   end  , ffend)
      flex  = c(fiflex,  FLX , ffflex)
      n_win = length(start)
      all_wins = 1:n_win
    }
  }



  if(verbose){
    digits = nchar(trunc(abs(len)))

    cat("\r\n Number of windows:",n_win)
    cat("\r\n Number of samples used:",(n_win-1)*inc+win)
    cat("\r\n Number of seconds used:",((n_win-1)*inc+win)/SR)
    cat0("\r\n Proportion of len used: ",((n_win-1)*inc+win)/len*100,"%\r\n")

    actual = 0;iterations = 0;
    while ( actual+win<=len) {
      iterations = iterations+1
      cat0('\r\nW',lead0(iterations, nchar(n_win))," | Samples: ", lead0(actual+1,digits),' - ',lead0(actual+win,digits),
           " | Seconds: ", lead0(actual+1/SR,digits),' - ',lead0((actual+win)/SR,digits))
      actual = actual+inc;
      #cat(actual+win,"\n\r")
    }
  }
  res = data.frame(window=all_wins, start_s=start,end_s=end,
                   mid_s=mid, duration = end-start+1, flex = flex,
                   start_t = tstart + (start-1)/SR, end_t = tstart + end/SR,
                   mid_t = tstart + (mid - 1)/SR )
  if(return=="number")  return(n_win)
  else if (return =="all") return(res)
}


#' Interpolate by window
#' Interpolate windowed data to original samplerate. 
#' useful to overlay computed indexes on original timeseries
#'
#' @param windowsList
#' @param winSec 
#' @param incSec 
#' @param SR 
#'
#' @return
#' @export
#'
#' @examples
winInter = function(windowsList, winSec, incSec, SR){
  warning("This function is a mess, please ask some developer to refactor it!")
  if(class(windowsList)!="list") {windowsList = list(windowsList)}
  #cat("Interpolating",incSec*SR,"sample between each HRV datapoint (linear) \r\n")
  inc=incSec*SR
  win= winSec*SR
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



#' Applies a function over moving windows
#'
#' @name byWin
#'
#' @param x a rats object or another object for which as.rats methods are available
#' @param WIN integer. the window size in the time unit of x (e.g. seconds)
#' @param INC integer. the increment size in the time unit of x (e.g. seconds)
#' @param flex logical. Creates additional shorter windows at the start and end
#' of the signal to maximize the coverage of the original series
#' @param position Should the time value of each calculation be positioned as the
#' first, middle or last sample of each window? 
#' @param FUN The function to be executed
#' @param ... Additional parameters to be passed to FUN
#' @param cores integer or "max". The number of cores to be used for parallelization.
#' if set to 1, parallelization is turned off. Very simple tasks can be faster if run
#' as single core due to the overhead of scheduling the task and returning the
#' result can be greater than the time to execute the task itself, resulting in
#' poor performance.
#'
#' @return 
#' @export
#'
#' @examples
byWin = function(x, WIN, INC, flex = TRUE, FUN, ..., 
                 position = c("mid", "first", "last"), cores = "max"){
  # DEBUG
  # stop("debug byWin")
  # x = lr10$all_CC_1$SC$s1
  # WIN = 300
  # INC= 0.1
  # FUN = "myRes"
  # flex=T
  # position = "mid"
  # l = list(newMin = 0.5, newMax = 0.8)
  # cores = 1
  #######
  x = as.rats(x)
  tu = if(is.rats(x)) timeUnit(x) else "unit of time"
  SR = frequency(x)
  win = WIN*SR
  win2 = round(win/2)
  inc = INC*SR
  position = match.arg(position)
  winz  = nwin(x, WIN, INC, SR, return = "all", flex = flex)
  n_win = nrow(winz)
  pos = switch (position,
    mid = winz$mid_s,
    start = winz$start_s,
    end = winz$end_s
  )
  l=list()
  l = list(...)
  
  if(cores == 1){
    cat(paste0("\r\nApplying ",
               deparse(substitute(FUN))," by ",WIN," ",tu, " windows, ",INC," ",tu," increments, in single core mode.\r\n")) #verified!
    FUN = match.fun(FUN)
    resList = vector(mode = "list", length = n_win)

    for(k in 1:n_win){
      ab = winz$start_s[k]:winz$end_s[k]
      resList[[k]] = do.call(FUN, c(list(x[ab]), l))
    }
    
  } else{
    ####### SETUP PARALLELIZATION
    # progresbar
    pb <- progress::progress_bar$new(
      format = "Calculation::percent [:bar] :elapsed | ETA: :eta",
      total = n_win,    #number of iterations
      width = 60, 
      show_after=0 #show immediately
    )
    
    progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar
    
    # allowing progress bar to be used in foreach -----------------------------
    progress <- function(n){
      pb$tick(tokens = list(letter = progress_letter[n]))
    } 
    
    opts <- list(progress = progress)
    
    
    #parallelization
    warningsFile = "MyWarnings"
    if (file.exists(warningsFile)) {
      unlink(warningsFile)
    }
    if(cores == "max")  cores=parallel::detectCores()-1
    
    # cat(paste0("\r\nPerforming parallelized computation of ",
    #            deparse(substitute(FUN)), " using ",cores," cores.\r\n")) 
    cat(paste0("\r\nApplying ",
               deparse(substitute(FUN))," by ",WIN," ",tu, " windows, ",INC," ",tu," increments, and using ",cores," cores.\r\n")) #verified!
    FUN = match.fun(FUN)
    cl <- parallel::makeCluster(cores[1], outfile=warningsFile) #not to overload your computer
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    pb$tick(0)
    resList <- foreach::foreach(
      k = 1:n_win,
      .options.snow = opts, .errorhandling='pass'
    ) %dopar% {
      ################# START OF THE PARALLEL JOB
      ab = winz$start_s[k]:winz$end_s[k]
      res = do.call(FUN, c(list(x[ab]), l))
      # res = FUN(x[ab], ...)
      
      ################# END OF THE PARALLEL JOB
    }
    parallel::stopCluster(cl) 
  }
  newTS = unlist(resList)
  
  
  res = rats(newTS, start = start(x)+round(pos[1]/SR), frequency = 1/INC, windowed = list(WIN, INC, flex = flex, winz))
  
  return(res)
}
#' @rdname byWin
wapply <- byWin


#' Approx by window
#' when you have a value for each window (even overlapping ones) and you
#' want to transform the values to a time series with the original sampling rate
#'
#' @param a a numerical vector, resulting from a windowing
#' @param winSec the size of the window used to calculate a, in seconds
#' @param incSec the increment of each window used to calculate a, in seconds
#' @param SR the number of samples per second of the destination time-series
#' @param midPoints if true, the value are considered to be at the center of the window, and the result is padded right by half window
#'
#' @return
#' @export
#'
#' @examples

# approxByWin = function(a, winSec, incSec, SR=frequency(x), midPoints=T){
#   # a = res
#   # winSec=4
#   # incSec = 1
#   # SR  = 100
#   warning("I think this function has issues!")
# 
#   inc=incSec*SR
#   win= winSec*SR
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


#' Rescale every sample of a time series over a moving window 
#' 
#' @description  This function allows to visualize the short-term dynamics of a
#' time series using the full vertical resolution of a screen
#' @param x numeric representation of a time-series. Ideally a rats object.
#' @param WIN numeric. Size of the window in the original series time unit (e.g. seconds)
#' @param newMin 
#' @param newMax 
#' @param cores 
#'
#' @export
rescaleByWin = function(x, WIN, newMin, newMax, cores=1){
  pointRescale = function(x){
    #x è un vettore, ma a me interessa solo la posizione centrale
    if(length(x)<3) stop ("poinRescale need longer objects")
    mid = round(length(x)/2)+1
    if(mid<2){ stop("bad mid point in pointrescale")}
    if(!is.numeric(x)) stop("in pointrescale x must be numeric")
    x = as.numeric(x) #strip classesand attributes
    ran = range(x)
    xx = x[mid]
    return((xx - ran[1])/(ran[2]- ran[1]))
    # return((newMax-newMin)/(ran[2]-ran[1]) * (xx - ran[2]) + newMax)
  }
  res = byWin(x, WIN=WIN, INC=1/frequency(x),flex = TRUE, FUN=pointRescale, cores=cores)
  res = res * (newMax-newMin) + newMin
}


#' Set argument/attribute to a list without overwriting existing values
#' 
#' This is most useful when dealing with ... and setting default behaviours that 
#' can be overweritten by ...
#'
#' @param arg 
#' @param value 
#' @param argList 
#'
#' @return

setArg = function(arg, value, argList){
  if(is.null(argList[[arg]])) argList[[arg]] = value
  argList
}


