# This new (old) approach, instead of getting a summary(e.g. median) for each occurrence of an epoch (or category),
# first pastes together the series correspoding to each epoch and THEN calculates the summarizing function
# no results on IM
#This approach may actually flatten up the results as given synchrony high variability,
# the median of a long segment tends to be the same, independently on how much synchrony was observed.
# a better summarizing metric with this approach could be % of time synchrony is above a given threshold



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
#' @param category 
#' @param categoryIndex 
#' @param artefact.rm 
#' @param mergeEpochs 
#' @param shift_start integer in seconds (e.g. -5 shifts all epoch starts by 5 seconds before )
#' @param shift_end integer in seconds (e.g. 7 shifts all epoch end by 7 seconds later )
#'
#' @return
#' @export
#'
#' @examples
epochSeries = function(x, signal, sync, series, category, categoryIndex,
                       mergeEpochs=FALSE, artefact.rm=TRUE, shift_start = 0, shift_end = 0){
  UseMethod("epochSeries", x)
}

#' @export
epochSeries.DyadExperiment = function(x, signal, sync, series, category, categoryIndex,
                                      mergeEpochs=FALSE, artefact.rm=TRUE, shift_start = 0, shift_end = 0){
  for(nSession in seq_along(x)){
    session = x[[nSession]]
    prog(nSession, length(x))
    x[[nSession]] = epochSeries(x[[nSession]], signal=signal, sync=sync, series=series,
                category=category, categoryIndex=categoryIndex, mergeEpochs=mergeEpochs,
                artefact.rm=artefact.rm, shift_start=shift_start, shift_end=shift_end )
  }

  return(x)
}

#' @export
epochSeries.DyadSession = function(x, signal, sync, series, category, categoryIndex, mergeEpochs, artefact.rm, shift_start, shift_end){
  # print(match.call())
  # print(missing(sync))
  # if(!missing(sync)) print(str(sync))
  ##DEBUG
  # x = d2[[1]]
  # signal="SC"
  # x$SC$artefacts = rbind(data.frame(start=c(160,300), end=c(200,400)),x[[signal]]$artefacts)
  # sync="amico2"
  # ziobilly = x$SC$amico2$sync
  # attributes(ziobilly) = NULL
  # x$SC$amico2$sync = rats(ziobilly, start=start(x$SC$amico2$sync)[1], frequency = frequency(x$SC$amico2$sync))
  # series = "sync"
  # category="IM"
  # categoryIndex="micro"
  # mergeEpochs = F
  # artefact.rm=T
  # shift_start = -2
  # shift_end = 2
  ###
  goodSyncs = names(x[[signal]])[sapply(x[[signal]],is.sync)]
  if(missing(sync)){

  }else if(! sync %in% names(x[[signal]])){
      stop ('sync was ',sync,' and should be one of: ',paste(goodSyncs,collapse = " "))
  }
  resName = paste0(c(category,"_",categoryIndex,"_",c(if(!missing(sync)){sync},series)),collapse = "")

  #select series according to sync and serieskey
  if(!missing(sync)){
    goodSeries = names(x[[signal]][[sync]])[sapply(x[[signal]][[sync]],is.rats)]
    if(length(goodSeries) == 0) stop ("'sync' argument must point to a list containing 'series'")
    if(!series %in% goodSeries) stop ('series was ',series,' and should be one of:',paste(goodSeries,collapse = " "))
    xseries = x[[signal]][[sync]][[series]]
  } else {
    goodSeries = names(x[[signal]])[sapply(x[[signal]],is.rats)]
    if(!series %in% goodSeries) stop ("series was ", series," and should be one of:",goodSeries )
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
  cate = x[[category]]
  #applica gli shift
  cate$start = cate$start + shift_start
  cate$end = cate$end + shift_end
  # print(str(cate))
  if(!is.factor(cate[[categoryIndex]])) stop("categoryIndex must be a factor column in the ",category,"'s epochs table. This are the values:", paste(colnames(cate[which(is.factor(cate))]),collapse = " ") )
  resList = list()
  #istanzia i contenitori vuoti per ciascun livello
  for(lev in levels(cate[[categoryIndex]])){
    # dur = sum((cate[cate[[categoryIndex]]==lev,"end"] - cate[cate[[categoryIndex]]==lev,"start"] )*frequency(xseries))
    if(mergeEpochs)
      resList[[lev]]=numeric()
    else
      resList[[lev]]=list()
  }
  names(resList) = levels(cate[[categoryIndex]])

  remSeries = xseries
  i=1
  for(i in 1:nrow(cate)){ #for each epoch
    if(!is.na(cate[[categoryIndex]][i]) && !is.null(cate[[categoryIndex]][i])) {
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


        #aggiungi la finestra al vettore di resList corrispondente al livello di categoryIndex
        win = window(xseries, start = cate$start[i], end = cate$end[i])
        win = c(padding, win)

        resList[[cate[[categoryIndex]][i]]] = c(resList[[cate[[categoryIndex]][i]]], list(win))
        names(resList[[cate[[categoryIndex]][i]]])[length(resList[[cate[[categoryIndex]][i]]])] = paste0(dyadId(x),sessionId(x),"|",cate$start[i], "-",cate$end[i])

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
#' @param experiment 
#' @param signal 
#' @param sync 
#' @param series 
#' @param category the category used to 
#' @param by currently unused. In future will be used to split the data by experimental group, participant, or any other relevant condition
#' @param FUN 
#' @param ... 
#'
#' @description  This function extracts the selected epochs from every session of a "DyadExperiment" object and puts them in
#'  a simple list. Categories must be created with epochSeries() beforehand.
#' @export
extractEpochs = function(experiment, signal, sync, series, category, categoryIndex, by, ...){
  if(!missing("by")) stop("by is not implemented yet.")
  UseMethod("extractEpochs",experiment) 
}

#' @export
extractEpochs.DyadExperiment = function(experiment, signal, sync, series, category, categoryIndex, by,   ...){
  if(missing(category) | missing(categoryIndex) | missing(series)){
      stop("category, categoryIndex, series, must all be specified")
  }
  epochsName = paste0(c(category,"_",categoryIndex,"_",c(if(!missing(sync)){sync},series)),collapse = "")
  resName =    paste0(c(category,"_",categoryIndex,"_",c(if(!missing(sync)){sync},series)),collapse = "")
  

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



