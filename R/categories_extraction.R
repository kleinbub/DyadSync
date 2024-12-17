# This new (old) approach, instead of getting a summary(e.g. median) for each occurrence of an epoch (or category),
# first pastes together the series correspoding to each epoch and THEN calculates the summarizing function
# no results on IM
#This approach may actually flatten up the results as given synchrony high variability,
# the median of a long segment tends to be the same, independently on how much synchrony was observed.
# a better summarizing metric with this approach could be % of time synchrony is above a given threshold


#' Join together all sessions of a given Diary
#'
#' @param exp a DyadExperiment object
#' @param diaryName string. The name of the Diary to be exported
#'
#' @return a dataframe containing all sessions of the selected Diary
#' @export
#'
#' @examples
oneDiary =function(exp, diaryName){
  if(!is.DyadExperiment(exp)) stop("x must be a DyadExperiment object")
  if(!diaryName %in% names(exp[[1]])) stop(paste("diaryName must be one of", paste(names(exp[[1]]),collapse = ", ") ))
  allb = c()
  for(i in 1:length(exp)){
    x = exp[[i]][[diaryName]]
    x$UID = UID(exp[[i]])
    allb = rbind(allb, x)
  }
  return(allb)
}


#' Prints a summary of all epochs
#' 
#' @param x a DyadExperiment or DyadSession object.
#' @param diaryName string. The name of the Diary to be summarized
#'
#' @return nothing.
#' @export
#'
#' @examples
catSummary = function(x, diaryName) {
  if(is.DyadExperiment(x)){
    x = oneDiary(x,diaryName)
  }else if(is.DyadSession(x)){
    x = x[[diaryName]]
  } else stop("x must be a DyadExperiment or DyadSession object.")
  s = summary(x)
  print(s)
  return()
}


#' Extract series chunks corresponding to epochs
#' 
#' epochSeries cycles through all sessions in an experiment.
#' For each session extracts a series, then cuts it
#' according to categories epochs
#'
#' @param x a DyadExperiment object
#' @param signal  string. The name of a DyadSignal present in x
#' @param sync If missing, series is searched in the signal (s1, s2, ...). If specified, series is searched within a sync object (PMBdev, CCFBest, ...). 
#' @param series 
#' @param diary 
#' @param category 
#' @param artefact.rm 
#' @param mergeEpochs 
#' @param shift_start numeric. Seconds to shift all start values (e.g. -5 shifts all epoch starts by 5 seconds before )
#' @param shift_end numeric. Seconds to shift all end values (e.g. 7.5 shifts all epoch end by 7.5 seconds later )
#' @param summaryFUN 
#' @param target one of "epoch", "before", "after". Determines whether to extract the actual epochs or only what precedes or follows it.
#' If it's not set to "epoch", shift_start or shift_end must be specified.
#' @param ... 
#' 
#' @return
#' @export
#'
#' @examples
epochSeries = function(x, signal, sync, series, diary, category,
                       mergeEpochs=FALSE, artefact.rm=TRUE, shift_start = 0, shift_end = 0, target=c("epoch", "before", "after"),
                       summaryFUN="median", ...){
  UseMethod("epochSeries", x)
}

#' @export
epochSeries.DyadExperiment = function(x, signal, sync, series, diary, category,
                                      mergeEpochs=FALSE, artefact.rm=TRUE, shift_start = 0, shift_end = 0,target=c("epoch", "before", "after"),
                                      summaryFUN="median", ...){
  for(nSession in seq_along(x)){
    session = x[[nSession]]
    prog(nSession, length(x))
    x[[nSession]] = epochSeries(x[[nSession]], signal=signal, sync=sync, series=series,
                diary=diary, category=category, mergeEpochs=mergeEpochs,
                artefact.rm=artefact.rm, shift_start=shift_start, shift_end=shift_end, target=target,
                summaryFUN = summaryFUN, ...)
  }

  return(x)
}

#' @export
epochSeries.DyadSession = function(x, signal, sync, series, diary, category,
                                   mergeEpochs, artefact.rm, shift_start, shift_end,target,
                                   summaryFUN, ...){
  # print(match.call())
  # print(missing(sync))
  # if(!missing(sync)) print(str(sync))
  #DEBUG
  # rm(list=ls())
  # library(DyadSync)
  # load("C:/Users/Kleinbub/OneDriveLink/__Ricerche/2021_biofeedback validation/BIOFEEDBACK LAB/SD_engine_v2.RData")
  # x = d2[[1]]
  # signal="SC"
  # sync="amico1"
  # series = "sync"
  # diary="SELF"
  # category="SELF"
  # mergeEpochs = F
  # artefact.rm=T
  # shift_start = 0
  # shift_end = 0
  # summaryFUN = "median"
  # target ="epoch"
  ###
  summaryFUN = match.fun(summaryFUN)
  target = match.arg(target, choices = c("epoch", "before", "after"))
  goodSyncs = names(x[[signal]])[sapply(x[[signal]],is.sync)]
  if(missing(sync)){

  }else if(! sync %in% names(x[[signal]])){
      stop ('sync was ',sync,' and should be one of: ',paste(goodSyncs,collapse = " "))
  }
  resName = paste0(c(diary,"_",category,"_",c(if(!missing(sync)){sync},"_",series)),collapse = "")

  #select series according to sync and serieskey
  if(!missing(sync)){
    goodSeries = names(x[[signal]][[sync]])[sapply(x[[signal]][[sync]],is.rats)]
    if(length(goodSeries) == 0) stop ("'sync' argument must point to a list containing 'series'")
    if(!series %in% goodSeries) stop ('series was ',series,' and should be one of:\r\n',paste(goodSeries,collapse = " "))
    xseries = x[[signal]][[sync]][[series]]
  } else {
    goodSeries = names(x[[signal]])[sapply(x[[signal]],is.rats)]
    if(!series %in% goodSeries) stop ("series was ", series," and should be one of:\r\n",goodSeries )
    xseries = x[[signal]][[series]]
  }

  ## remove Artefacts windows from xseries, if any
  if( artefact.rm && nrow(x[[signal]]$artefacts)>0 ){
    for(i in 1:nrow(x[[signal]]$artefacts)){
      # cat("\r\n",x[[signal]]$artefacts$start[i], " ", x[[signal]]$artefacts$end[i])
      window(xseries, start=x[[signal]]$artefacts$start[i],end=x[[signal]]$artefacts$end[i]) <- NA
    }

  }
  # seleziona il dataframe della categoria e crea un oggetto per ciascun livello
  cate = x[[diary]]
  
  
  if(target == "before"){
    #in questo caso l'end della finestra corrisponde allo start
    #e lo start corrisponde allo start + shift_start che deve essere negativo
    if(missing(shift_start) || shift_start >= 0) stop("With target='before', shift_start must be a non-zero negative number")
    if(!missing(shift_end) && shift_end != 0) message("Whit target='before' shift_end should be used with awareness!")
    cate$end   =  cate$start
  } else if(target == "after"){
    #in questo caso lo start della finestra corrisponde all'end
    #l'end corrisponde alla end + shift_end 
    if(missing(shift_end) || shift_end <= 0) stop("With target='after', shift_end must be a non-zero positive number")
    if(!missing(shift_start) && shift_start != 0) message("Whit target='after' shift_start should be used with awareness!")
    cate$start = cate$end
  }
  #applica gli shift in tutti e tre i casi
  cate$start = cate$start + shift_start
  cate$end = cate$end + shift_end
  cate$delta = cate$end - cate$start #ricalcola dopo le modifiche fatte
  # print(str(cate))
  if(!is.factor(cate[[category]])) stop("category must be a column of type 'factor' in the ",diary,"'s epochs table.\r\nThis are the possible values:\r\n", paste(colnames(cate[sapply(cate,is.factor)]),collapse = ", ") )
  resList = list()
  #istanzia i contenitori vuoti per ciascun livello
  for(lev in levels(cate[[category]])){
    # lev = levels(cate[[category]])[1]
    # dur = sum((cate[cate[[category]]==lev,"end"] - cate[cate[[category]]==lev,"start"] )*frequency(xseries))
    if(mergeEpochs)
      resList[[lev]]=numeric()
    else
      resList[[lev]]=list()
  }
  names(resList) = levels(cate[[category]])

  remSeries = xseries
  i=1
  for(i in 1:nrow(cate)){ #for each epoch
    if(!is.na(cate[[category]][i]) && !is.null(cate[[category]][i])) {
      if(cate$start[i]>=end(xseries)[1]){
        warning("In session ", dyadId(x),"-",sessionId(x), ", start of window ",i,": was equal to or beyond the series end.", call.=F)
        # lres[[i]] = NA
      } else { #if start is before the end of xseries, as it should...

        #if (by applying shift_start) cate$start is before the signal start, create a NA padding
        if(cate$start[i]<start(xseries)){
          padding = rats(start = cate$start[i],
                     end = start(xseries), 
                     frequency = frequency(xseries),
                     timeUnit=timeUnit(xseries), unit=unit(xseries))
          cate$start[i] = start(xseries)
        } else padding = NULL

        #if end goes beyond the xseries duration...
        if(cate$end[i] > end(xseries)[1] ){
          message("In session ", dyadId(x),"-",sessionId(x), ", end of window ",i,": was reduced to the series end.\n")
          cate$end[i]= end(xseries)[1]
        }

        # rimuovi la finestra dallo series di segnale residuo remSeries:
        window(remSeries, start = cate$start[i], end = cate$end[i]) <- NA ## broken?


        #aggiungi la finestra al vettore di resList corrispondente al livello di category
        win = window(xseries, start = cate$start[i], end = cate$end[i])
        if(!is.null(padding)){
          win = c(padding, win)
        }
        cn = paste0(c(sync, series),collapse="_")
        x[[diary]][i,cn] = summaryFUN(as.numeric(na.omit(win)), ...)

        resList[[cate[[category]][i]]] = c(resList[[cate[[category]][i]]], list(win))
        names(resList[[cate[[category]][i]]])[length(resList[[cate[[category]][i]]])] = paste0(dyadId(x),sessionId(x),"|",cate$start[i], "-",cate$end[i])

      }
    }
  }
  if(mergeEpochs){
    resList[["remaining"]] = list(na.omit(remSeries))
    trueResList = resList
    resList = trueResList
    for(i in 1:length(resList)){
      if(length(resList[[i]])>0)
        resList[[i]] = do.call(c, resList[[i]])
    }
  } else {
    remSeries = as.numeric(remSeries)
    resList[["remaining"]] = list(remSeries)
}
  #save object
  x[[signal]][[resName]] = resList
  x
}



#' @title extract epochs in a simple list
#'
#' @param experiment 
#' @param signal 
#' @param sync 
#' @param series 
#' @param diary the Diary 
#' @param category the Category of the Diary to be extracted
#' @param by currently unused. In future will be used to split the data by experimental group, participant, or any other relevant condition
#' @param ... currently unused.
#' 
#' @description  This function extracts the selected epochs from every session of a "DyadExperiment" object and puts them in
#'  a simple list. Categories must be created with epochSeries() beforehand.
#' @export
extractEpochs = function(experiment, signal, sync, series, diary, category, by, ...){
  if(!missing("by")) stop("by is not implemented yet.")
  UseMethod("extractEpochs",experiment) 
}

#' @export
extractEpochs.DyadExperiment = function(experiment, signal, sync, series, diary, category, by,   ...){
  if(missing(diary) | missing(category) | missing(series)){
      stop("diary, category, series, must all be specified")
  }
  epochsName = paste0(c(diary,"_",category,"_",c(if(!missing(sync)){sync},"_",series)),collapse = "")

  #check names
  keepNames = unique(unlist(lapply (experiment, function(session){
    goodNames = names(session[[signal]])[!sapply(session[[signal]],is.sync) & !sapply(session[[signal]],is.rats)]
    if(! epochsName %in% goodNames) stop(epochsName, " was not found in session. Have you run epochSeries() beforehand? Do you need to specify sync? Found names: ", paste0(goodNames,collapse=" ") )
    
    names(session[[signal]][[epochsName]])
  })))
  keepNames = sort(keepNames)
  #controlla se mergeEpochs era T o F
  mergeEpochs = !any(unlist(lapply (experiment, function(session){
    lapply(session[[signal]][[epochsName]],is.list)
  })))
  #instanzia una lista vuota con keepNames elementi
  res = vector("list", length(keepNames))
  names(res) = keepNames

  for(session in experiment){
    if(!is.null(session[[signal]][[epochsName]])){
      sesNames= names(session[[signal]][[epochsName]])
      for(sesName in sesNames){
        newElement = session[[signal]][[epochsName]][[sesName]]
        newElement = if(mergeEpochs  || is.list(newElement)) newElement else list(newElement)
        res[[sesName]] = c(res[[sesName]], newElement)
      }
    }
  }
  # boxplot(res)
  # 
  # unlist(lapply(res,FUN, ...))
  res
}



