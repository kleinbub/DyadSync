##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################ 
## Credits
# Author: Johann R. Kleinbub
# Contact: johann.kleinbub@gmail.com
############################################################################################################ 

##Data importer juggernaut

############################################################################################################
############################################################################################################
############################################################################################################

#' Title
#'
#' @param path 
#' @param signalNames a vector of string defining the name of the signals (same order as s1Col, s2Col)
#' @param start an optional vector of integers, or strings in mm:ss format which specify the time (in seconds) of the first observation of each file.
#' Note that by default time series will be set with start = 0, which is different from the default value of \code{\link[stats]{ts}}.
#' @param end an optional vector of integers, or strings in mm:ss format which specify the time (in seconds) of the last observation of each file.
#' If the data is longer or shorter, remove tail
#'  or pads with zeroes. Useful to cut signals when a session is over or to equalize lengths. The default value,
#'  FALSE, keeps the original length of the signal
#' @param duration an optional vector of integers, or strings in mm:ss format which specify the time (in seconds) of the duration of each file.
#' If the data is longer or shorter, remove tail or pads with zeroes. Only one between end or duration should be specified.
#' @param pairBind if true, every two consequent files in path get matched to a single DyadSession, useful if the data of participant
#'                 1 and 2 are stored in separate files.
#' @param winTerpolate a list containing two numeric values, winSec and incSec. If the data to be read has been generated through a moving
#'                     window filter, it can be reversed to the given SR by setting the size (in seconds) of the window and its increment.
#' @param namefilt a string used to select only specific files in a given path
#' @param idOrder either NA or a character vector. idOrder is used to interpret the filenames to correctly identify and label the cases. 
#'   If NA, no attempt to identify cases will be done. The character vector may contain one or more of the following strings in any given order: "id" (unique identifier of a dyad),
#'   "session" (session number if the same dyad has multiple sessions),"group" (if the data must be devided across groups or experimental conditions),
#'   "role" (if pairBind = TRUE, to assign the role of different files), or "x" (a placeholder to skip irrelevant information in filenames)
#'   The strings can be abbreviated.
#' @param idSep character vector (or object which can be coerced to such) containing regular expression(s).
#'   If idOrder is not NA, this will be used as separator to split the filenameas and identify "id", "session", and "group"
#'   informations.
#' @param ... additional options to be passed to read.table (es: skip, header, etc)
#' @param s1Col 
#' @param s2Col 
#' @param s1Name 
#' @param s2Name 
#' @param unit character. A descriptor of the unit of measurement of the series values.
#' @param SR 
#' @param timeUnit character. A descriptor of the time unit of the series values.
#' Defaults to seconds as in "samples per second"
#'
#' @return
#' @export
#'
#' @examples
readDyadSignals = function(
                      path,  #location and format of experiment files
                      s1Col,s2Col, #one or multiple columns identifiyng the patient and clinician for each signal
                      s1Name, s2Name,
                      signalNames, #how to rename each signal (same order as s1Col, s2Col)
                      SR,    #sampling rate at which the data is acquired
                      start, #
                      end,  #remove tail or pads with zeroes. Useful to cut signals when a session is over or to equalize lengths
                      duration,
                      unit, timeUnit="seconds", 
                      pairBind = F, #if true, each two files in the path get matched to a single DyadSession, useful if the data of patient
                                    #and clinician are saved on separate files
                      winTerpolate = list(winSec=NULL,incSec=NULL), #if data comes from a moving windows analysis it should be
                                                                    #restored to the original data sampling rate
                      namefilt = NA,
                      idOrder= c(), # c("id","session","group","role"), #the order of identifiers in the filenames. role is only used with pairbind = T
                      idSep = "_", #the separator of identifiers in the filename
                      ... #additional options to be passed to read.table (es: skip, header, etc)
                      ){
####debug #####
  # path = data_d
  # maxSeconds=F
  # s1Col=c(2)
  # s2Col=c(2)
  # signalNames = c("SC")
  # ###signalNames = c("start", "stop", "LFHF", "PNN50", "HF", "RRMean",  "RRsd")
  # ###winTerpolate = list(winSec = 30, incSec = 5)
  # winTerpolate = list(winSec = NULL, incSec = NULL)
  # 
  # 
  # namefilt = NA
  # pairBind=T
  # idOrder= c("id","session","x","r")
  # pairBind=F
  # idOrder= c("id","session","x")
  # idSep = "_"
  # options=list("skip"=0, "sep"="\t", "header"=F)
  # 
  # s1Col=c(3)
  # s2Col=c(4)
  # s1Name = "patient"
  # s2Name = "therapist"
  # signalNames = c("SC")
  # start = c(vsCC)
  # duration = c(durCC)
  # SR=100
  # idOrder = c("s","x","x","x","id")
  # unit = "uS"
  # timeUnit = "seconds"
######
  if(!missing(end) && !missing(duration)) stop("only one between end and duration must be specified")
  if(missing(start)) start = 0
  if(length(s1Col)!=length(s2Col) || length(s1Col)!=length(signalNames)) stop ("s1Col, s2Col, signalNames must have same length")
  
  # imp = genericIO(path,namefilt,idOrder,idSep, pairBind)
  imp = genericIO(path,namefilt,idOrder,idSep, pairBind, ...)
  lf= imp$lf
  sess= imp$sess
  dyadIds = imp$dyadIds
  group = imp$group
  role = imp$role
  filenames = imp$filenames
  shortNames = imp$shortNames
  ndyads = imp$ndyads
  nFiles = imp$nFiles
  
  
  if(!is.numeric(unlist(lf[[1]][s1Col]))) stop("s1Col column is not numeric. Maybe there is a header? Set skiprow to 1 or more?")
  if(!is.numeric(unlist(lf[[1]][s2Col]))) stop("s2Col column is not numeric. Maybe there is a header? Set skiprow to 1 or more?")


  len <- lapply(lf, function(x) length(x[[1]]))
  #check if some file are shorter than 50% of the mean length and remove them 
  removeFile = unlist(lapply(seq_along(len), function(i){if(len[[i]] / mean(unlist(len)) <0.5) {
      warning("File ",i,': ',shortNames[i]," is shorter (",timeMaster(round(len[[i]]/SR),out="m"),
              "s) than the 50% of the mean length", call. = F);return(i)}
    } ))
  # if(length(removeFile)>0){
  #   lf[removeFile] = len[removeFile] = dyadIds[removeFile] = group[removeFile] = sess[removeFile] = NULL
  #   filenames =filenames[-removeFile]; shortNames = shortNames[-removeFile]
  #   
  # }



  if(pairBind){
    if(!isTRUE(all.equal(s1Col,s2Col))) {warning("If pairBind is true, you most probably want to have s1Col and s2Col to be equal! Watch what your're doing!", call. = F)
    }  else {
      s2Col = s2Col + ncol(lf[[1]])
    }
    
    cat0("\r\nCombining '",levels(factor(unlist(role)))[1],"' and '", levels(factor(unlist(role)))[2],"' files\r\n")
    #add 'role' labels to each column. 'lr' = list renamed
    lf =  Map(function (x,i) {names(x)=paste(names(x),role[[i]],sep="_"); x}, lf, seq_along(lf) )
    
    #combine dyad's members in a single dataframe. 'lrb' = list renamed bound
    lfu = lf[seq(1,length(lf)-1, by=2)]
    lfu = Map(function(i,j){
      if(dyadIds[[i]] != dyadIds[[i+1]]) stop ("file ",i,":",dyadIds[[i]], "is different from file",i+1,": ", dyadIds[[i+1]])
      if(role[[i]] == role[[i+1]]) stop ("role of file ",i,": ",role[[i]]," is equal to file ",i,": ",role[[i+1]])
      cat0("merging ",dyadIds[[i]],"_",sess[[i]]," ", role[[i]]," (as patient) with ",dyadIds[[i+1]],"_",sess[[i]]," ", role[[i+1]]," (as clinician)\r\n" )
      lfu[[j]] = cbind(lf[[i]],lf[[i+1]])
    },seq(1,length(lf)-1, by=2), 1:(length(lf)/2) )
    lf = lfu
    rm(lfu)
    names(lf) = unlist(dyadIds[c(TRUE,FALSE)])
    dyadIds = dyadIds[c(TRUE,FALSE)]
    group = group[c(TRUE,FALSE)]
    sess = sess[c(TRUE,FALSE)]
  }
  
  #check if any file has only NA's
  removeFile = unlist(lapply(seq_along(lf), function(i){
    if(any(apply(lf[[i]][c(s1Col, s2Col)], 2, function(k) all(is.na(k)) ))) {
      warning("In dyad ",dyadIds[[i]],' - session ',sess[[i]],", at least one column was all NA's. The whole session was removed", call. = F)
      return(i)
    }

  } ))
  if(length(removeFile)>0){
    lf[removeFile] = len[removeFile] = dyadIds[removeFile] = group[removeFile] = sess[removeFile] = NULL
    filenames =filenames[-removeFile]; shortNames = shortNames[-removeFile]
    
  }
  
  
  if(!is.null(winTerpolate$winSec) & !is.null(winTerpolate$incSec)){
    cat("\r\nInterpolating files\r\n")
    prog(1,2)
    lf=winInter(lf,winTerpolate$winSec,winTerpolate$incSec,SR)
    prog(2,2)
  }  
  
  #start checks
  start = timeMaster(start,out="sec")
  if(any(start==1)) warning("start = 1 is usually wrong. Use default value of zero to read a series from the beginning.")
  if(length(start) == 1 && ndyads>1){
    message("start = ",start," was used for all ",ndyads," dyads")
    start=rep(start,ndyads)
  } else if(length(start) != ndyads){
    stop("start must be a single value or be provided for each ",ndyads ,"file. Only ",length(start)," values provided:\r\n",start)
  }
  
  #end checks
  if(!missing(end)){ #if end is NOT missing
    end = timeMaster(end,out="sec")
    duration = end - start
  }  

  #if duration is NOT missing, or was set through 'end' procede to trim files 
  if(!missing(duration) || !missing(end)){ 
    duration = timeMaster(duration, out="sec")
    cat("\r\nTrimming files (samples)\r\n")
    
    if(length(duration) == 1 && ndyads>1) {
      duration=rep(duration,ndyads)
      message("duration was used for all ",ndyads," dyads")
    } else if(length(duration) !=ndyads) stop("end/duration must be provided for each file",length(duration) ,"!=",ndyads)
    
    #check if any duration is NA and replace it with the actual file duration
    for(i in seq_along(lf)){
      if(is.na(duration[i])) duration[i] = nrow(lf[[i]])
    }

    
    print(data.frame("file"=shortNames, "original"=sapply(lf,nrow),"destination"=duration*SR,"NAs added"=sapply(seq_along(lf), function(i){
      if(nrow(lf[[i]])<duration[i]*SR) "*" else "-"
    })),row.names = F)
    
    #check if some file are  shorter than duration and add NAs
    lf <- lapply(seq_along(lf), function(i){
      if(nrow(lf[[i]])<duration[i]*SR){
        lf[[i]][(nrow(lf[[i]])+1):(duration[i]*SR),] = NA
      }
      return(lf[[i]]) })
    #resize file according to settings
    lf <- lapply(seq_along(lf),function(i){lf[[i]][1:(duration[i]*SR) ,]})
  }
  ndyads = length(lf)
  
  ###############################
  ## reading report
  len <- lapply(lf, function(x) length(x[[1]]))
  outtable =data.frame(
    "dyad" = unlist(dyadIds),
    "session" = unlist(sess),
    "group" = unlist(group),
    "start" = timeMaster(round(start),out="min"),
    "duration" = timeMaster(floor(unlist(len)/SR), out="min"),
    #"max length" = timeMaster(end, out="min"),
    # "Filename" = unlist(shortNames),
    
    row.names = if(!pairBind) shortNames else shortNames[c(T,F)]
  )
  # if(!pairBind) outtable$role = NULL
  # if(length(end)==1 && !end) outtable$max.length = NULL
  print(outtable)
  cat("\r\n",ndyads, "dyads successfully imported from",nFiles,"files.\r\n")
  ################################
  
  
  #Populates the objects
  sessions = vector(mode="list",length=length(lf))
  iSession = 1
  for(iSession in 1:length(lf)){
    session = lf[[iSession]]
    signalList = lapply(seq_along(s1Col), function(i) {
      DyadSignal(name=signalNames[i],
                 s1=rats(session[,s1Col[i]],frequency=SR,
                         start=start[iSession], timeUnit=timeUnit, unit=unit),
                 s2=rats(session[,s2Col[i]],frequency=SR,
                         start=start[iSession], timeUnit=timeUnit, unit=unit),
                 SR = SR, s1Name = s1Name, s2Name = s2Name,
                 sessionId=sess[[iSession]],
                 dyadId=dyadIds[[iSession]],
                 groupId=group[[iSession]])
    })
    #generates sessions
    ses = DyadSession(groupId = group[[iSession]], sessionId= sess[[iSession]], dyadId = dyadIds[[iSession]],
                      signalList=signalList, s1Name = s1Name, s2Name = s2Name, fileName = shortNames[[iSession]] )
    sessions[[iSession]] = ses
  }
  experiment = DyadExperiment(path, sessions)
  
  # experiment = DyadExperiment(path,
  #   Map(function(session,nSession){
  #   #for each type of signal, add a new DyadSignal to the present DyadSession
  #   #These are defined as s1Col paTer pairs.
  #   signalList = lapply(seq_along(s1Col), function(i) {
  #     DyadSignal(name=signalNames[i],
  #                s1=rats(session[,s1Col[i]],frequency=SR,
  #                        start=start[nSession], timeUnit="second", unit=unit),
  #                s2=rats(session[,s2Col[i]],frequency=SR,
  #                        start=start[nSession], timeUnit="second", unit=unit),
  #                SR = SR, s1Name = s1Name, s2Name = s2Name,
  #                sessionId=sess[[nSession]],
  #                dyadId=dyadIds[[nSession]],
  #                groupId=group[[nSession]])
  #   })
  #   #generates sessions
  #   ses = DyadSession(groupId = group[[nSession]], sessionId= sess[[nSession]], dyadId = dyadIds[[nSession]],
  #                     signalList=signalList, s1Name = s1Name, s2Name = s2Name, fileName = shortNames[[nSession]] )
  #   return(ses)
  #   
  # },lf,seq_along(lf))
  # )
#   if(pairBind){
#     names(experiment$sessions) = names(lf)
#   } else {
    names(experiment) = paste(group,dyadIds,sess,sep="_")
#   }
  
  #lr=lapply(lr, na.omit)
  # cat("\r\n\r\nReport:\r\n")
  # print( data.frame("names"=names(lf),
  #                   "orig_duration_min"=timeMaster(as.numeric(unlist(len))/SR, out="min"),
  #                   "final_duration_min"= timeMaster(as.numeric(unlist(lapply(lf, function(x) length(x[,s2Col[1]]))))/SR, out="min"),
  #                   "orig_samp_size"=as.numeric(unlist(len)),
  #                   "final_samp_size"= as.numeric(unlist(lapply(lf, function(x) length((x[,s2Col[1]])))))
  # )
  # )
  # cat("\r\n\r\nInitial skipped rows are not considered in this report.\r\n")
  return(experiment)
  
}






