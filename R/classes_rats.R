#'###############################################################################
#'                                _              _                      
#'                      _ __ __ _| |_ ___    ___| | __ _ ___ ___ 
#'                     | '__/ _` | __/ __|  / __| |/ _` / __/ __|
#'            /\//\//\/  | | (_| | |_\__ \ | (__| | (_| \__ \__ \
#'           |/\//\//\/|_|  \__,_|\__|___/  \___|_|\__,_|___/___/
#'                                                               
#'###############################################################################
#'             _   _               _   _   _                         _        
#'    _ _ __ _| |_(_)___ _ _  __ _| | | |_(_)_ __  ___   ___ ___ _ _(_)___ ___
#'   | '_/ _` |  _| / _ \ ' \/ _` | | |  _| | '  \/ -_) (_-</ -_) '_| / -_|_-<
#'   |_| \__,_|\__|_\___/_||_\__,_|_|  \__|_|_|_|_\___| /__/\___|_| |_\___/__/
#'                                                                            
#'###############################################################################
#'
#' State of the file:
#' this is a v0.2 let's say, in which instead of having a $x $y sructure, I am using
#' a y and attr(y, 'time') format. this is good because x * 10 is preserved. 

#' An ideal time series object to represent pyhisio signals should have the
#' following information:
#' start_date, 
#' start, end, duration (milliseconds)
#' sampling per second (Hz)
#' windowing history: size, inc, flex, windows
#' PRINCIPLES:
#' * rats are legion: there is no singular form "rat" as there is no singular "timeserie" # nolint
#' * rats are tiny: 'rats' is never capitalized.
#' * rats are social: rats(start=0, end=10, f=1) and rats(start=10, end=20, f=1) do
#'   not overlap and can be combined to a single series starting at 0 and ending at 20.
#'   Time intervals are half-open intervals, including the start value
#'   but excluding the end one.
#' * rats are familiar: syntax is intuitive for ts() or zoo() users
#' * Holes fit in cheese, but not in rats: a rats' length must correspond to the rats'
#'   duration multiplied by the frequency. Missing data in rats must
#'   be represented with NAs.
#'   

roundFast = function(x) {
  as.integer(x + sign(x)*0.5)
}

#' @title ~rats: Rational Time Series
#' @description rats is the creator for a nimble S3 class of regular time series.
#' rats are similar to \link[stats]{ts} or \link[zoo]{zooreg} classes, but these
#' classes have quirks related to quarterly periods and other econometric
#' analytic traditions, while rats is more oriented toward engineering and computer
#' science applications and implements the logics of these disciplines.
#' 
#' Features:
#' \itemize{
#'    \item fractional frequency or period is perfectly supported
#'    \item time intervals are defined as half-open intervals [A, B)
#'    \item robust and non-destructive subsetting/extraction of series using the time coordinates with functions such
#'    as \link[DyadSync]{window.rats}, \link[DyadSync]{'window<-.rats'}, \link[DyadSync]{at}
#'    \item robust and non-destructive manipulation, expansion, and combination of series with functions such as
#'    \link[DyadSync]{c.rats} and \link[DyadSync]{'window<-.rats'}
#'    \item practical tools to explore(\link[DyadSync]{'print.rats'}, \link[DyadSync]{'str.rats'}) and
#'    plot (\link[DyadSync]{'plot.rats'}) rats
#' }
#' 
#' 
#' @param data vector. A vector of any type representing the data of the series
#' @param start numeric. The beginning of the series in time, in cycle units.
#' @param end numeric. The end of the series in time, in cycle units.
#' @param duration numeric. Duration of the series in cycle units.
#' @param frequency numeric. Either frequency of period must be present to specify the
#' sampling rate of the signal. Frequency is espressed in number of observations
#' per cycle. While period is the duration of each observation in fractions of
#' cycle units. The period is equal to 1/frequency.
#' @param period numeric. 
#' @param windowed list of three elements winSec, incSec, flex, describing eventual
#' windowing procedures used to create the series. 
#' @param timeUnit character. A descriptor of the frequency cycle unit. Used only
#' in print functions.
#' @param unit character. A descriptor of the unit of measurement of the series values. Used only
#' in print functions.
#' 
#' @details 
#' In the ~rats class, time intervals are defined as half-open intervals
#' with start and end and window methods including the start value but not
#' the end value
#' This is also referred to as [A, B)
#' or "inclusive start, exclusive end". 
#' See: \href{https://www.cs.utexas.edu/users/EWD/transcriptions/EWD08xx/EWD831.html}{Dijkstra, E.W. (1982). Why numbering should start at zero}
#' 
#' 
#' Thus, given a series starting at 0, with 5000 samples of information at 1000 samples/second
#' you only have samples for time 0, 1, 2, ... until 4999 so this would be represented as
#' 00:00 - 00:04
#' time(1) == 0  == 00:00
#' time(N) == 4999 == 00:04-1*period
#' 
#' In another example, rats(start=0, end=10, f=1) and rats(start=10, end=20, f=1) do not overlap.
#' 
#' The class works by storing the y values as a vector, and the x, or "time" values as an attribute,
#' together with other relevant metadata
#'   
#'
#' @return a rats time series, which is a vector of any type with
#' @export
#'
#' @examples
#' 
#' 
#' 
rats = function(data, start=0, end, duration, frequency=1, period,
                windowed = list(size=NULL, increment=NULL, flex=NULL, table=NULL),
                timeUnit="cycle", unit=NULL) {
  # print(match.call())
  # ######debug
  # data = 7.6
  # start=1024
  # # end=1
  # frequency = 10
  # # duration = end-start
  # stop("debug rats")
  # #############
  
  
  fd = FALSE
  if(missing(period)){
    #flag to remember if frequency was just estimated (lower priority)
    if(missing(frequency)){fd = T }
    #frequency c'è di default
    period = 1/frequency  
  } else {
    if(missing(frequency)){
      frequency = 1/period
    } else {
      #both are specified. Are they coherent?
      if(period != 1/frequency){
        stop("Period must be equal to 1/frequency.")
      }
    }
  }
  
  #set durations:
  if(!missing(duration)){
    if( missing(end))                    end = start + duration
    if( missing(start) && !missing(end)) start = end - duration
    if(duration != end - start) stop("~rats error code 1: Specified arguments led to incoherent series duration")
    if((!missing(data) && length(data)>0) && !fd ){
      alt_d = signif(length(data)/frequency,6)
      if(duration != alt_d) stop("~rats error code 2: Specified arguments led to incoherent series duration")
    }
    
  } else {
    #if 'duration' is missing
    durations = c()
    if(!missing(end) && !missing(start)) durations = c(durations, round(end-start,digits = 10))
    if(!missing(data) && length(data)>0) durations = c(durations, round(length(data)/frequency, digits = 10))
    if(length(unique(durations))>1)  stop("~rats error code 3: Specified arguments led to incoherent series duration")   
    duration = durations[1]
    if( missing(end)) end = start + duration
    if( missing(start) && !missing(end)) start = end - duration
    
  }
  
  #finally if you have a good duration and the data you can get frequency
  if(!missing(data) && length(data)>0 && fd){
    frequency =  length(data)/duration
  }
  
  
  
  
  #generate the time values
  if(length(start)  >0 && !is.na(start)  &&
     length(end)    >0 && !is.na(end)    &&
     length(period) >0 && !is.na(period)
  ){
    #' this old strategy to generate x values was prone to bugs
    # delta = end-period
    # #fix a bug the sign of by is wrong
    # if(delta < start){
    #   if(isTRUE(all.equal(start, delta))){
    #     delta = start
    #   } else {
    #     stop("error in rats. Line 190. sequence created negative numbers")
    #   }
    # }
    # x = seq(start,delta,by=period)
    
    #' whith this strategy you should have the guarantee o having exactly 1 x value
    #' for each y one, no floating point errors, and other seq() quirks.
    #' The first value is always going to be == start, and the following ones are
    #' increments of 1 period.
    #' This approach is also 10x times faster!
    nrep = roundFast((end - start)/period-1)
    x = start + cumsum(c(0,rep(period, times= nrep)))

  } else{
    x = numeric()
  }
  n = length(x)
  
  if(!missing(data) && length(data)>0){
    # if(duration*frequency != length(data)) stop("duration must equal data/frequency")
    if(length(data)>n) data= data[1:n]
    if(length(data)<n) stop("Not enough data was provided for the specified duration.")
    y = data
  } else {
    #generate empty data
    y = as.numeric(rep(NA, n))
  }
  res = y
  
  if(!missing(windowed) && !all(sapply(windowed, is.null))){
    if(!is.list(windowed)) stop("windowed must be a list")
    if(length(windowed) <3) stop("windowed must contain at least:
                                 window size, window increment, flex")
    if(length(windowed) >4) stop("too many elements provided in the windowed argument")
    WIN = windowed[[1]]
    INC = windowed[[2]]
    FLX = windowed[[3]]
    if(length(WIN)>1 || !is.numeric(WIN) || is.null(WIN)) stop("window size must be numeric of length 1")
    if(length(INC)>1 || !is.numeric(INC) || is.null(INC)) stop("window increment must be numeric of length 1")
    if(length(FLX)>1 || !is.logical(FLX) || is.na(FLX))   stop("window flex must be logical of length 1")
    if(length(windowed) == 4){
      TAB = windowed[[4]]
    } else {
      TAB = nwin(y, WIN, INC, frequency, return = "all", flex = FLX)
    }
    windowed = list("size"=WIN, "increment"=INC, "flex"=FLX, table= TAB)
  }
  
  attributes(res) = c(attributes(res),
                      list("x"=round(x, digits = 10),
                           "start" = start,
                           "end" = end,
                           "duration" = duration,
                           "frequency" = frequency,
                           "period" = period,
                           "windowed" = as.list(windowed),
                           "n" = n,
                           "timeUnit" = timeUnit,
                           "unit" = unit
                      ))
  class(res) = "rats"
  return(res)
}

#' @export
is.rats = function(x){inherits(x,"rats") && length(x)}

#' Subsets rats by time
#' Intervals always include the 'start' and exclude the 'end' values.
#' Windows can alternatively be specified through start, end, or one among
#' start and end and a duration.
#' 
#'
#' @param x 
#' @param start numeric. Start of the window
#' @param end numeric. End of the window
#' @param duration numeric. The duration of the window.
#'
#' @return a subset rats object
#' @export
#'
#' @examples
#' x = rats(1:100, start=0, frequency = 10)
#' 
#' 
#' window(x, start=2, duration=2)
#' 
#' #Replace values within a window
#' window(x, start=0.21, end=0.5) = 99 
#' 
#' 
# win = window(xseries, start = cate$start[i], end = cate$end[i])
# x = xseries
# start = cate$start[i]
# end = cate$end[i]
window.rats = function(x, start, end, duration){
  
  #se tutti e tre, devono essere coerenti
  if(!missing(start) && !missing(end) && !missing(duration)){
    if(end-start != duration) stop("Duration must be equal to end - start.")
  }
  if(!missing(start) && !missing(end)){
    duration = end - start
  } else if(!missing(start) && !missing(duration)){
    end = start + duration
  } else if(!missing(end) && !missing(duration)){
    start = end - duration
  } else if(!missing(duration)){
    stop("Duration must always be specified together with start or end.")
  } else {
    #using attr here as a weird hack because if start [argument] is missing
    #then start [function] is not being called :-0
    if (missing(start)) start = attr(x,"start")
    if (missing(end))   end   = attr(x,"end")
    duration = end - start
  }
  if(length(start(x))>1) stop("'rats' were improperly declassed to 'ts'. Try restarting the session. Or contact devs if this repeats.")
  
  if(round(duration, digits=10) < round(period(x), digits=10)){
    duration = period(x)
    end = start + period(x)
    warning("zero length duration was coerced to 1 sample")
  }

  # due to floating point errors even simple operations (e.g. subtraction)
  # with non-integer numbers can lead to unpredictable errors 
  # eg:  1024.1 - 0.1 < 1024  [TRUE!]
  # print(1024.1 - 0.1, digits = 18) [1023.99999999999989]
  # with rounding this is fixed.
  # print(round(1024.1 - 0.1,digits = 10), digits = 18)
  
  start = round(start, digits=10)
  end   = round(end, digits=10)
  duration = round(duration, digits=10)
  x_time = x$x#round(x$x, digits = 10)
  y = x$y

  # start parameter can be whatever fractional number
  #Actual start instead must be a discrete number of periods before or after the
  # original signal start.
  n_delta_pre = trunc(round((start(x)-start)/period(x), digits=8))
  new_start = start(x) - n_delta_pre*period(x)
  
  n_delta_post = floor(round((end - end(x))/period(x), digits=8))
  new_end = end(x) + n_delta_post*period(x)
  
  #this should NEVER happen
  if(new_end - new_start< .Machine$double.eps^0.5){
    warning ("windowing returned unexpected zero length series. Please contact developers")
    return(NULL)
  }

  
  #generate the new time scale from new_start to new_end
  #end is non inclusive, so the last x value is new_end-1*period
  nrep  = roundFast((new_end - new_start)/period(x)-1)
  new_x  = new_start + cumsum(c(0,rep(period(x), times = nrep)))
  new_x  = round(new_x, digits = 10)
  new_y  = rep(NA, length(new_x))

  
  #if the original signal is inside the new range, write the old values.
  keep_i = which(x_time %in% new_x)
  if(length(keep_i)>0){
    where_i = which(new_x %in% x_time)
    if(!all.equal(new_x[where_i], x_time[keep_i])) warning("window.rats encountered logical error 1.")
    new_y[where_i] = y[keep_i]
  }
  res = rats(new_y, start = new_start, end = new_end, frequency = frequency(x),
             windowed = attr(x, "windowed"), timeUnit = timeUnit(x), unit = unit(x)
  )
  if(!isTRUE(all.equal(new_x, res$x))) warning("window.rats encountered logical error 2.")
  
  return(res)

  
}

#' @export
"window<-.rats" = function(x, start, end, duration, values){
  
  ##DEBUG
  # x = lr10$all_CC_1$SC$s1
  # start = 120
  # end = 220.1
  # a = window(x, start, end)
  # values = 1:1001
  # stop("debug")
  #############
  
  #se tutti e tre, devono essere coerenti
  if(!missing(start) && !missing(end) && !missing(duration)){
    if(end-start != duration) stop("Duration must be equal to end - start.")
  }
  if(!missing(start) && !missing(end)){
    duration = end - start
  } else if(!missing(start) && !missing(duration)){
    end = start + duration
  } else if(!missing(end) && !missing(duration)){
    start = duration - end
  } else if(!missing(duration)){
    stop("Duration must always be specified together with start or end.")
  } else {
    if (missing(start)) start= start(x)
    if (missing(end))   end  = end(x)
    duration = end - start
  }
  
  if(round(duration, digits=10) <period(x)){
    duration = period(x)
    end = end+period(x)
    warning("zero length duration was coerced to 1 sample")
  }
  
  start = round(start, digits=10)
  end   = round(end, digits=10)
  duration = round(duration, digits=10)
  x_time = round(x$x, digits = 10)
  y = x$y
  
  # start parameter can be whatever fractional number
  #Actual start instead must be a discrete number of periods before or after the
  # original signal start.
  n_delta_pre = trunc(round((start(x)-start)/period(x), digits=10))
  new_start = start(x) - n_delta_pre*period(x)
  
  n_delta_post = floor(round((end - end(x))/period(x), digits=10))
  new_end = end(x) + n_delta_post*period(x)
  
  #generate the new time scale from new_start to new_end
  #end is non inclusive, so the last x value is new_end-1*period
  new_x  = new_start + cumsum(c(0,rep(period(x), round((new_end - new_start)/period(x)-1,digits=10))))
  new_x  = round(new_x, digits = 10)
  
  #check provided values
  if(length(values)==1) values = rep(values, length(new_x))
  if(length(values) !=length(new_x)) stop(paste(length(values), "values provided, but",length(new_x),"where expected from the window range."))
  
  #in assigment mode, all the real data must be kept
  fin_start = min(start(x), new_start)
  fin_end   = max(end(x),   new_end)
  fin_x  = fin_start + cumsum(c(0,rep(period(x), round((fin_end - fin_start)/period(x)-1,digits=10))))
  fin_x  = round(fin_x, digits = 5)
  fin_y  = rep(NA,length(fin_x)) 
  
  
  #write the old values inside the new range.
  where_i = which(fin_x %in% x_time)
  fin_y[where_i] = y
  
  #write the new values inside the new range.
  where_i = which(fin_x %in% new_x)
  fin_y[where_i] = values
  
  res = rats(fin_y, start = fin_start, end = fin_end, frequency = frequency(x),
             windowed = attr(x, "windowed"), timeUnit = timeUnit(x), unit = unit(x)
  )
  if(!all.equal(fin_x, res$x)) warning("window.rats encountered logical error 3.")
  
  return(res)
}

#' Get rats value at given time
#'
#' @param x a rats object
#' @param time a time in a format accepted by \link[DyadSync]{timeMaster}
#'
#' @return a rats object of length 1
#' @export
#'
#' @examples
#' at(rats(500:100, frequency=10), time = 12.378)
at = function(x, time){
  if(!is.rats(x))stop("at requires a ~rats")
  if(!is.numeric(time)) time = timeMaster(time, out="sec", digits=7)
  window(x, start=time, duration=period(x))
}

#' Get rats sample number at given time
#'
#' @param x a rats object
#' @param time a time in a format accepted by \link[DyadSync]{timeMaster}
#' @details this can be useful e.g. after windowing
#' @return the index of the vector x where x$x matches time
#' @export
#'
#' @examples
#' at_s(rats(500:100, frequency=10), time = 12.378)
at_s = function(x, time){
  if(!is.rats(x))stop("at requires a ~rats")
  if(!is.numeric(time)) time = timeMaster(time, out="sec", digits=7)
  match = window(x, start=time, duration=period(x))
  which(abs(x$x - match$x) < 1e-6  )
}

#' @export
"[.rats" = function(x,i){
  if(all(is.na(i))) return(NA_real_)
  if(all(is.logical(i))){
    #all good
  } else if(all(is.numeric(i))){
    if(max(i)>attr(x,"n")) stop("Subset out of bounds")
    if(min(i)<1) stop("0 or negative subscripts are not supported")
  } else {stop("Subscripting must be numeric or logical")}

  if(length(i)>attr(x,"n")) stop("Subset was longer than the data. Use window.rats for flexible expansion.")
  
  if(length(i)==0) return(numeric(0))
  x1 = .subset(time(x), i) #time, or x
  x2 = .subset(x, i)       #values, or y
  
  
  wa = attr(x, "windowed")
  if(!is.null(wa$table)){
    nt = wa$table[wa$table$mid_t >= x1[1] &
                    wa$table$mid_t <= x1[length(x1)], ]
    wa$table = nt
  }
  
  res = rats(x2,start=x1[1],frequency = frequency(x), 
             timeUnit = timeUnit(x), unit = unit(x),
             windowed = wa)
  return(res)
  
}

#' @export
"[<-.rats" = function(x,i,values){
  if(all(is.na(i))) return(x)
  if(length(i)>attr(x,"n")) stop("Subset out of bounds")
  if(length(i)==0) {
    return(x)
  } else {
    k = unclass(x)
    k[i] = values
    class(k) = "rats"
    return(k)
  }
  
}
#' @export
"$.rats" = function(x,i){
  i = match.arg(i, c("x","time","index","y","values"))
  if(i %in% c("x","time","index"))
    return(time(x))
  else if(i %in% c("y","values")){
    attributes(x) = NULL
    return(x)
  } else
    stop("Admissible values for rats$ are: $x, $time, $index for time(rats), or $y, $values for the series values as vector")
}

#' @export
length.rats = function(x){
  if(attr(x,"n") != length(unclass(x))) stop("Actual length and n attribute mismatch")
  attr(x,"n")
}

#' @export
"$<-.rats" = function(x,i,values){
  stop("Direct edit of ~rats are not allowed. Use window() or create new rats")
}

#' @export
str.rats = function(x, vals = 5, digits=4, ...){
  cat0("~~rats ")
  if(length(attr(x,"start")) >0 && !is.na(attr(x,"start")))
    cat0("from ",attr(x,"start"), " to ",attr(x,"end"))
  else cat0("of length 0")
  if(attr(x,"frequency") >=1){
    cat0(" | ",attr(x,"frequency")," samples/",attr(x, "timeUnit"))
  } else {
    if(attr(x,"period") == 1)
      cat0(" | 1 sample every ", attr(x, "timeUnit"))
    else
      cat0(" | 1 sample every ",attr(x,"period")," ", attr(x, "timeUnit"),"s")
    
  }
  if(!is.null(attr(x,"windowed")[[1]])){
    winprint = attr(x,"windowed")[1:3]
    cat0(" | windowed: ",paste(winprint,collapse=" "))
  }
  cat0(": ")
  #print them all
  if(length(time(x))==0){
    print("#@BUG: this is not implemented for some reason. Please notify the authors.")
  } else {
    realn = min(vals,length(x))
    toprint = 1:realn
    valz = x[toprint]
    if(is.numeric(valz)) valz = round(valz, digits)
    # k = rbind("values"=valz,"time"=round(time(x)[toprint],digits))
    cat(valz)
    if(length(x)>vals) cat("...")
  }
  cat0("\r\n")
}

#' Print a rats object
#' 
#' @param x 
#' @param vals numeric. How many values should be printed?
#' @param digits numeric. How many digits should be used to approximate values
#'
#' @export
print.rats = function(x, vals=20, digits= 4){
  cat0("~~rats time series ")
  cat0("[1:",length(x$y),"] ")
  if(length(attr(x,"start")) >0 && !is.na(attr(x,"start")))
    cat0("from ",attr(x,"start"), " to ",attr(x,"end"))
  else cat0("of length 0")
  if(attr(x,"frequency") >=1){
    cat0(" | ",attr(x,"frequency")," samples/",attr(x, "timeUnit"))
  } else {
    if(attr(x,"period") == 1)
      cat0(" | 1 sample every ", attr(x, "timeUnit"))
    else
      cat0(" | 1 sample every ",attr(x,"period")," ", attr(x, "timeUnit"),"s")
    
  }
  if(!is.null(attr(x,"windowed")[[1]])){
    winprint = attr(x,"windowed")[1:3]
    cat0(" | windowed: ",paste(winprint,collapse=" "))
  }
  cat0(".")
  if(is.numeric(vals) && length(time(x)) > vals){
    #number of shown values must be even and greater than 2
    if(vals < 2) vals = 2
    if(vals %% 2 != 0) vals = vals - 1
    
    #print the first half of the values
    timeDecimals = -log(attr(x,"period"),10)
    toprint = 1:(vals/2)
    k1 = rbind("values"=x[toprint]$y,"time"=round(time(x)[toprint], digits))
    if(is.numeric(k1)) k1[1,] = as.character(round(k1[1,], digits))
    k1[2,] = signif(as.numeric(k1[2,]), )
    
    toprint = (length(time(x))-vals/2+1):length(time(x))
    k2 = rbind("values"=x[toprint],"time"=round(time(x)[toprint], digits))
    if(is.numeric(k2)) k2[1,] = as.character(round(k2[1,], digits))
    
    hid = length(time(x)) - vals
    if(hid>1000)hid = paste0(floor(hid/1000),"K")
    k = cbind(k1,paste0("\U2026[",hid," more]\U2026"),k2)
    
  } else {
    #print them all
    if(length(time(x))==0){
      toprint = numeric()
      k = rbind("values"="numeric(0)","time"="numeric(0)")
      
    } else {
      toprint = 1:length(time(x))
      valz = x[toprint]
      if(is.numeric(valz)) valz = round(valz, digits)
      k = rbind("values"=valz,"time"=round(time(x)[toprint],digits))
    }
    
  }
  k = as.data.frame(k)
  colnames(k)=NULL
  print(k)
}


#' @export
start.rats = function(x){attr(x,"start")}
    
#' @export
end.rats   = function(x){attr(x,"end")}
#' @export
duration = function(x){UseMethod("duration",x)}
#' @export
duration.rats = function(x){attr(x,"duration")}
#' @export
frequency.rats = function(x){attr(x,"frequency")}
#' @export
period = deltat.rats = function(x){attr(x,"period")}
#' @export
cycle.rats = function(x){(time(x)%%frequency(x))%%1*frequency(x)+1}
#' @export
time.rats = function(x){attr(x,"x")}
#' @export
timeUnit = function(x){UseMethod("timeUnit",x)}
#' @export
timeUnit.rats = function(x){attr(x,"timeUnit")}
#' @export
unit = function(x){UseMethod("unit",x)}
#' @export
unit.rats = function(x){attr(x,"unit")}

#' @export
windowed = function(x){UseMethod("windowed",x)}
#' @export
windowed.rats = function(x){
  wa = attr(x,"windowed")
  res = wa
  if(!is.null(wa[[1]])){
    res = wa[[4]]
    attributes(res) = c(attributes(res),
                        list("windowed"=wa[1:3]),
                        timeUnit = timeUnit(x)
    )
    class(res) = c("winTable",class(res))
  }
  return(res)
}

#' @export
print.winTable = function(x, n=10L, ...){
  wa = attr(x,"windowed")
  tu = attr(x, "timeUnit")
  if(is.null(wa[[1]])){
    cat0("This rats is not the output of a windowing function")
  } else {
    cat0("This rats was generated through a windowing procedure using window
sizes of ", wa[[1]], " ",tu,
         "s and increments of ", wa[[2]], " ",tu, "s")
    if(wa[[3]]) cat0(" using flex windows")
    cat0(".\r\nTable of windows:\r\n")
    N = min(n,nrow(x))
    TAB = as.data.frame(x[1:N,])
    print(TAB, ...)
    if(nrow(x)>N) cat("...")
  }
}

#' Combine rats
#' 
#' Various rats objects can be combined if their x values don't overlap, and
#' if they share the same frequency and timeUnit.
#'
#' @param ... 
#'
#' @return a rats object.
#' @export
#'
#' @examples
#' 
#' a1 = rats(runif(10),start=20, frequency=10)
#' a2 = rats(letters[1:22], start=0, end=2.1, frequency = 10)
#' a3 = rats(rep(0,13),start = 21, frequency = 10, timeUnit = "second")
#' 
#' l = list(a1,a2,a3)
#' print(l)
c.rats = function(...){
  l = list(...)
  ## DEBUG
  # l = list(rats(1:10, start=3.5), rats(1:10, start = 20))

  if(length(l)==1) return(l[[1]])
  
  #if there are some non rats, ratsify them
  nr = which(sapply(l,\(x){!is.rats(x)}))
  rr = which(sapply(l,\(x){is.rats(x)}))
  if(length(nr)>0) {
    for(i in nr){
      if(class(l[[i]]) %in% c("ts","zooreg")){
        l[[i]] = as.rats(l[[i]])
      } else {
        #find the first rat behind it
        previous = rr[which(rr<i)]
        if(length(previous) > 0){
          example = l[[previous[length(previous)]]]
          old = l[[i]]
          attributes(old) =NULL
          l[[i]] = rats(old, start = end(example), frequency = frequency(example),
                        timeUnit = timeUnit(example), unit = unit(example),
                        windowed = attr(example, "windowed")
          )
        } else {
          #there were no previous rats
          stop("combining non rats with rats in this way is not implemented yet")
        }
      }
    }
  }
  lfreq = unique(sapply(l,frequency))
  if(length(lfreq)>1) stop("All rats must have the same frequency.")
  lunit = unique(sapply(l,timeUnit))
  if(length(lunit)>1) stop("All rats must have the same cycle unit.")
  lvunit = unique(sapply(l,unit))
  if(length(lvunit)>1) stop("All rats must have the same unit of measurement.")
  lwind = unique(lapply(l,function(x){attr(x, "windowed")[1:3]}))
  if(length(lwind)>1) stop("All rats must have the same windowed attribute.")
  lwind = lwind[[1]]
  
  ends   = unlist(lapply(l,end))
  starts = unlist(lapply(l,start))
  #Sort the rats
  ord = order(starts)
  ends = ends[ord]
  starts = starts[ord]
  l = l[ord]
  #check if all rats are contiguous
  deltas =  starts[2:length(starts)] - ends[1:(length(ends)-1)]
  #negative deltas mean overlap. Bad rats!
  round(deltas, digits=10)
  if(any(deltas<0))stop("Rats' heads and tails must not overlap.")
  #positive deltas mean gaps. Create empty rats to fill them. Poor, empty, rats...
  if(any(deltas>0)){
    #one rat for each gap
    for(i in which(deltas>0)){
      l = c(l, list(rats(start=ends[i], end = starts[i+1],
                         frequency=frequency(l[[1]]))) )
    }
    #Recursion magic. Praise the rat-gods.
    l = do.call("c",l)
    
  } else {
    #Finally join all rats' y together into a final... y.
    finaly = c(unlist(lapply(l,function(a) a)))
    #The rat is now whole!
    return(rats(finaly, start=starts[1], frequency = lfreq,
                windowed = lwind,timeUnit = lunit, unit = lvunit))
  }
  
}

#' @export
# @BUG range su rats da come minimo 1 non so perché
range.rats = function(..., na.rm=FALSE){
  c(min(..., na.rm=na.rm), max(...,na.rm=na.rm))
  
}

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
na.omit.rats = function(x){
  #if nas are at beginning or end return a rats, else return cheese
  k = unclass(x)
  
  nna = which(!is.na(x))
  tk = time(x)[nna]
  k = k[nna]
  
  if(length(unique(diff(tk)))==1) {
    #nas where only at the beginning and end
    k = rats(k, start = tk[1], frequency = frequency(x),
             windowed = attr(x,"windowed"), timeUnit = timeUnit(x),
             unit = unit(x) )
  } else {
    attributes(k) = list(
      "x" = tk,
      "na.removed" = which(is.na(x)),
      start=tk[1],
      unit = unit(x)
    )
    class(k) = "cheese"
  }
  
  return(k)
}

#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.rats = function(x, ...){
  l  = list(...)
  
  if(is.null(l[["type"]])) l[["type"]]="l"
  if(is.null(l[["xlab"]])) l[["xlab"]]=paste0("Time (",timeUnit(x),")")
  if(is.null(l[["ylab"]])) l[["ylab"]]=paste0("Values (",unit(x),")")
  
  k = list(x=x$x,y=x$y)
  l = c(k,l)
  if(all(is.na(k$y))) k$y = as.numeric(x$y)
  if(!is.numeric(k$y)) stop("Only numerical rats can be plotted right now.")
  do.call("plot",l)
}

#' @export
lines.rats = function(x, ...){
  lines(x$x, x$y, ...)
}

#' @export
as.rats = function(x){
  UseMethod("as.rats",x)
}
#' @export
as.rats.rats = function(x){
  x
}
#' @export
as.rats.ts = function(x){
  xstart = start(x)
  if(length(xstart)==2){
    xstart = xstart[1]+ (xstart[2]-1)* 1/frequency(x)
  } else if(length(startx)==1){
    xstart = xstart[1]
  } else stop("unknown error")
  k = x
  attributes(k) = NULL
  rats(k, start=xstart, frequency = frequency(x))
}
#' @export
as.rats.zooreg = function(x){
  k = x
  attributes(k) = NULL
  rats(k, start=start(x), frequency = frequency(x))
}
#' @export
as.rats.numeric = function(x){
  k = x
  attributes(k) = NULL
  rats(k, start=0, frequency = 1)
}



#' @export
as.ts.rats = function(x){
  k = x
  attributes(k) = NULL
  ts(k, start=start(x), frequency = frequency(x))
}
#' @export
as.zooreg.rats = as.zoo.rats = function(x) {
  k = x
  attributes(k) = NULL
  zooreg(k,order.by =time(x), frequency = frequency(x))
}

#' @export
quantile.rats = function(x,...){quantile(x$y,...)}

#' @export
as.data.frame.rats = function(x,  ...) 
{ as.data.frame(x$y, ...)
}

#' @export
diff.rats = function(x, loseFrom=c("mid", "head","tail"), ...) {
  x1 = diff(as.numeric(x), ...)
  delta = (length(x) - length(x1))*period(x)
  loseFrom = match.arg(loseFrom, choices = c("mid", "head","tail"))
  newStart = if(loseFrom=="mid"){ start(x) + delta/2
  } else if(loseFrom=="tail") { start(x)
  } else {start(x) + delta}
  x1 = rats(x1, start=newStart, frequency = frequency(x))
  attr(x1, "windowed") = attr(x, "windowed")
  attr(x1, "timeUnit") = attr(x, "timeUnit")
  attr(x1, "unit") = attr(x, "unit")
  return(x1)
}
