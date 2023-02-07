################################################################################
#                                _              _                      
#                      _ __ __ _| |_ ___    ___| | __ _ ___ ___ 
#                     | '__/ _` | __/ __|  / __| |/ _` / __/ __|
#            /\//\//\/  | | (_| | |_\__ \ | (__| | (_| \__ \__ \
#           |/\//\//\/|_|  \__,_|\__|___/  \___|_|\__,_|___/___/
#                                                               
################################################################################
#             _   _               _   _   _                         _        
#    _ _ __ _| |_(_)___ _ _  __ _| | | |_(_)_ __  ___   ___ ___ _ _(_)___ ___
#   | '_/ _` |  _| / _ \ ' \/ _` | | |  _| | '  \/ -_) (_-</ -_) '_| / -_|_-<
#   |_| \__,_|\__|_\___/_||_\__,_|_|  \__|_|_|_|_\___| /__/\___|_| |_\___/__/
#                                                                            
################################################################################
#'
#' State of the file:
#' this is a v0.2 let's say, in which instead of having a $x $y sructure, I am using
#' a y and attr(y, 'time') format. this is good because x * 10 is preserved. 

#' An ideal time series object to represent pyhisio signals should have the
#' following information:
#' start_date, 
#' start, end, duration (milliseconds)
#' sampling per second (Hz)
#' windowing history: size, inc, flex
#' 
#' given a STS starting are at 0000, if you have 5000ms of information
#' you only have samples for time 0, 1, 2, ... until 4999 so this would be represented as
#' 00:00 - 00:04
#' Indeed time intervals are defined as half-open intervals
#' with start and end and window methods including the start value but not
#' the end value
#' 
#define the rats class rational time series

#' Rational Time Series
#' rats is the creator for a nimble S3 class of regular time series.
#' rats are similar to \link[stats]{ts} or \link[zoo]{zooreg} classes, but these
#' classes have many quirks related to quarterly periods and other econometric
#' analytic traditions, while rats is more oriented towards physical signals.
#' The crucial difference is that for rats, time intervals are half-open intervals,
#' which means that they include the start value but not the end one. This way
#' rats(start=0, end=10, f=1) and rats(start=10, end=20, f=1) do not overlap.
#' Secondarily, rats are stored as a list with $x and $y values, storing respectively
#' the temporal information, and the values. These latter can be numeric or not.
#' 
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
#' * rats are legion: there is no singular form "rat" as there is no singular "timeserie"
#' * rats are small: rats is never capitalized.
#' * rats fit well together: rats(start=0, end=10, f=1) and rats(start=10, end=20, f=1) do
#'   not overlap and are combined to a . Time intervals are half-open intervals, including the start value
#'   but excluding the end one.
#' * rats are familiar: syntax is intuitive for ts() or zoo() users
#' * Holes fit in cheese, but not in rats: a rats' length must correspond to the rats'
#'   duration multiplied by the frequency. In other terms missing data in rats must
#'   be represented with NA
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
rats = function(data, start=0, end, duration, frequency=1, period,
                windowed = list(winSec=NULL, incSec=NULL, flex=NULL),
                timeUnit="cycle", unit=NULL){
  
  # ######debug
  # start=1.9
  # end=10
  # frequency = 10
  # duration = end-start
  # #############
  if(missing(period)){
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

  durations = c()
  if(!missing(duration)) durations = c(durations,duration)
  if(!missing(start) && !missing(end)) durations = c(durations, end-start)
  if(!missing(data) && !missing(frequency) && length(data)>0){
    #qua puoi farlo se hai frequency
    durations = c(durations, signif(length(data)/frequency,6))}
  if(length(unique(durations))>1) stop("duration mismatch")   
  if(length(unique(durations))==0) stop("Insufficient information to build rats") 
  duration = unique(durations)
  
  if(!missing(start) &&  missing(end)) end = start + duration
  if( missing(start) && !missing(end)) start = end - duration


  #generate the time values
  x = seq(start,end-period,by=period)
  n = length(x)
  
  if(!missing(data) && length(data)>0){
    # if(duration*frequency != length(data)) stop("duration must equal data/frequency")
    if(length(data)>n) data= data[1:n]
    if(length(data)<n) stop("Not enough data was provided for the specified duration.")
    y = data
  } else {
    #generate empty data
    y = rep(NA, n)
  }
  res = y
  
  attributes(res) = c(attributes(res),
                      list("x"=x,
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
    start = duration - end
  } else if(!missing(duration)){
    stop("Duration must always be specified together with start or end.")
  } else {
    if (missing(start)) start= start(x)
    if (missing(end))   end  = end(x)
    duration = end - start
  }
  # stop ("at least two between start, end, and duration must be specified")
  
  if(signif(duration*frequency(x),6)%%1 != 0) stop("Duration must be a multiple of sampling rate")
  if(start < start(x) | end > end(x)) stop("The window was outside of the rats data boundary")
  
  cutter = which(time(x) >= start & time(x) < end )
  if(length(cutter)==0) stop("No data found in the given window")
  rats(x[cutter], start=start, frequency = frequency(x), windowed = attr(x, "windowed"),
       timeUnit = timeUnit(x), unit = unit(x))
  
}

#' @export
"window<-.rats" = function(x, start, end, duration, values){
  #se tutti e tre, devono essere coerenti
  if(!missing(start) && !missing(end) && !missing(duration)){
    if(end-start != duration) stop("Duration must equal to end - start.")
  }
  if(!missing(start) && !missing(end)){
    duration = end - start
  } else if(!missing(start) && !missing(duration)){
    end = start + duration
  } else if(!missing(end) && !missing(duration)){
    start = duration - end
  } else stop ("at least two between start, end, and duration must be specified")
  
  if(start < start(x) | end > end(x)) stop("The window was outside of the rats data boundary")
  
  cutter = which(time(x) >= start & time(x) < end )
  if(length(cutter)==0) stop("No data found in the given window")
  x[cutter] = values
  x
}

#' @export
"[.rats" = function(x,i){
  if(length(i)>attr(x,"n")) stop("Subset out of bounds")
  x1 = .subset(time(x), i)
  x2 = .subset(x, i)
  rats(x2,start=x1[1],frequency = frequency(x), windowed = attr(x, "windowed"),
       timeUnit = timeUnit(x), unit = unit(x))
}

#' @export
"[<-.rats" = function(x,i,values){
  if(length(i)>attr(x,"n")) stop("Subset out of bounds")
  k = unclass(x)
  k[i] = values
  class(k) = "rats"
  k
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
str.rats = function(x){
  str(unclass(x))
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
  cat0("from ",attr(x,"start"), " to ",attr(x,"end"))
  if(attr(x,"frequency") >=1){
    cat0(" | ",attr(x,"frequency")," samples per ",attr(x, "timeUnit"))
  } else {
    if(attr(x,"period") == 1)
      cat0(" | 1 sample every ", attr(x, "timeUnit"))
    else
      cat0(" | 1 sample every ",attr(x,"period")," ", attr(x, "timeUnit"),"s")
    
  }
  if(!is.null(attr(x,"windowed")[[1]])){
    cat0(" | windowed: ",paste(attr(x,"windowed"),collapse=" "))
  }
  cat0(".")
  if(is.numeric(vals) && length(time(x)) > vals){
    #number of shown values must be even and greater than 2
    if(vals < 2) vals = 2
    if(vals %% 2 != 0) vals = vals - 1
    
    #print the first half of the values
    timeDecimals = -log(attr(x,"period"),10)
    toprint = 1:(vals/2)
    k1 = rbind("values"=x[toprint],"time"=round(time(x)[toprint], digits))
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
    toprint = 1:length(time(x))
    valz = x[toprint]
    if(is.numeric(valz)) valz = round(valz, digits)
    k = rbind("values"=valz,"time"=round(time(x)[toprint],digits))
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
  if(length(l)==1) return(l[[1]])
  
  #if there are same non rats, ratify them
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
                        timeUnit = timeUnit(example)
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
  lvunit = unique(sapply(l,timeUnit))
  if(length(lvunit)>1) stop("All rats must have the same unit of measurement.")
  lwind = unique(lapply(l,function(x){attr(x, "windowed")}))
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
  print(l)
  if(all(is.na(k$y))) k$y = as.numeric(x$y)
  if(!is.numeric(k$y)) stop("Only numerical rats can be plotted right now.")
  do.call("plot",l)
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
    xstart = xstart[1]+ (xstart[2]-2)* 1/frequency(x)
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


# microbenchmark::microbenchmark(
#   ts  (1:100,start=0, frequency = 10)[1:10],
#   rats(1:100,start=0, frequency = 10)[1:10]
# )
# profvis::profvis({
#   x = rats(1:100,start=0, frequency = 10)
#   x[1:10] = NA
# 
# })

#' 
#' 
#' #examples
#' ##EXAMPLE
#' rm(list=ls())
#' x = rats(100:120, 0, freq=365, val="uS", tim="years")
#' x
#' print(x,vals = "all", digits=3)
#' str(x)
#' plot(x, col="red")
#' x
#' a2 = rats(letters[1:22], start=0, end=2.1, frequency = 10)
#' a1 = rats(sample(1:100,10),start=20, frequency=10) #end: 21
#' a3 = rats(rep("test",13),start = 21, frequency = 10, timeUnit = "second")
#' 
#' asd = c(a1,a2,a3)
#' print(asd, 50)
#' 
#' 
#' 
#' x = rats( start=0, duration=10, frequency = 0.5)
#' print(x)
#' str(x)
#' library(zoo)
#' 
#' #empty TS
#' x = rats( start=0, duration=10, frequency = 0.5)
#' y = zoo(order.by = time(x),frequency = 0.5 )
#' z = ts(start=0, end=10, frequency = 0.5)
#' length(x);length(y);length(z)
#' print(x);
#' print(y);
#' print(z)
#' 
#' na.omit(x)
#' na.omit(y)
#' na.omit(z)
#' 
#' x = c(x,rats(99,start=end(x), frequency = 0.5))
#' y = c(y,zoo(99,order.by =end(y)+1, frequency = 0.5))
#' z = c(z,ts(99,start=end(z), frequency = 0.5))
#' 
#' na.omit(x)
#' na.omit(y)
#' na.omit(z)
#' 
#' #full TS
#' #' frequency 0.5
#' #' y value: | a | b | c | d | e | f | g | h | i | j |
#' #' sample:  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10|
#' #' t value: 0 | 2 | 4 | 6 | 8 | 10| 12| 14| 16| 18| 20|
#' #'
#' rm(x,y,z)
#' data = letters[1:6]; freq=5
#' x = rats(data, start=0, frequency = freq)
#' z = ts  (data, start=0, frequency = freq)
#' y = zoo (data, order.by = seq(0,length(data)*freq, by=1/5),frequency = freq )
#' length(x);length(y);length(z)
#' print(x);
#' print(y);
#' print(z)
#' 
#' start(x);start(y);start(z)
#' end(x);end(y);end(z)
#' 
#' x = list(x, rats(letters[6:10], start=10, end=20, frequency = 0.5))
#' z = list(z, ts  (letters[6:10], start=10, end=20, frequency = 0.5))
#' y = list(y, zoo (letters[6:10], order.by = seq(10,20,by=2),frequency = 0.5 ))
#' length(x);length(y);length(z)
#' print(x);
#' print(y);
#' print(z)
#' 
#' start(x);start(y);start(z)
#' end(x);end(y);end(z)
#' 
#' 
#' 
#' # small frequency
#' #' a ts of 4 seconds, sampled 4 times per second
#' #' frequency 4
#' #' y value: | a | b | c | d | e | f | g | h | i | j |
#' #' sample:  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10|
#' #' t value: 0 | 2 | 4 | 6 | 8 | 10| 12| 14| 16| 18| 20|
#' #'
#' rm(x,y,z)
#' x = rats(letters[1:20], start=0, end=2, frequency = 10)
#' z = ts  (letters[1:20], start=0, end=2, frequency = 10)
#' y = zoo (letters[1:20], order.by = seq(0,10,by=1/10),frequency = 10 )
#' length(x);length(y);length(z)
#' str(x);
#' str(y);
#' str(z)
#' 
#' start(x);start(y);start(z)
#' end(x);end(y);end(z)
#' 
#' window(x, start=0.5, end=1.5)
#' time(window(y, start=0.5, end=1.5))
#' window(z, start=0.5, end=1.5)
#' 
#' 
#' # basic frequency
#' #' a ts of 4 seconds, sampled 4 times per second
#' #' frequency 4
#' #' y value: | a | b | c | d | e | f | g | h | i | j |
#' #' sample:  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10|
#' #' t value: 0 | 2 | 4 | 6 | 8 | 10| 12| 14| 16| 18| 20|
#' #'
#' rm(x,y,z)
#' x = rats(letters[1:20], start=0, end=20, frequency = 1)
#' z = ts  (letters[1:20], start=0, end=20, frequency = 1)
#' y = zoo (letters[1:20], order.by = seq(0,20,by=1),frequency = 1 )
#' length(x);length(y);length(z)
#' str(x);
#' str(y);
#' str(z)
#' 
#' start(x);start(y);start(z)
#' end(x);end(y);end(z)
#' 
