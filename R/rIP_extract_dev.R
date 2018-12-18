# This new (old) approach, instead of getting a summary(e.g. median) for each occurrence of an epoch (or category),
# first pastes together the stream correspoding to each epoch and THEN calculates the summarizing function
# no results on IM




#' Extract stream chunks corresponding to epochs
#' 
#' epochStream cycles through all sessions in an experiment.
#' For each session extracts a stream, then cuts it
#' according to categories epochs
#'
#' @param x 
#' @param signal 
#' @param sync 
#' @param streamKey 
#' @param category 
#' @param groupIndex 
#' @param artefact.rm 
#'
#' @return
#' @export
#'
#' @examples
epochStream = function(x, signal= "SC", sync="PMBest", streamKey="sync", category, groupIndex, mergeEpochs=F, artefact.rm="T"){
  UseMethod("epochStream", x)
}

#' @export
epochStream.DyadExperiment = function(x, signal, sync, streamKey, category, groupIndex, mergeEpochs, artefact.rm){
  res = Map(function(session, nSession){
    prog(nSession, length(x))
    print(nSession)
    epochStream(session, signal, sync, streamKey, category, groupIndex,mergeEpochs, artefact.rm )
  },x, seq_along(x))
  classAttr(res) = classAttr(x)
  res
}

#' @export
epochStream.DyadSession = function(x, signal, sync, streamKey, category, groupIndex, mergeEpochs, artefact.rm){
  if(!sync %in% names(x[[signal]])) stop ('sync should be one of:',names(x[[signal]]))
  resName = paste0(c(toupper(category),"_",totitle(c(sync,streamKey))),collapse = "")
  
  #select stream according to sync and streamkey
  if(sync!="none"){
    if(!streamKey %in% names(x[[signal]][[sync]])) stop ('sync should be one of:',names(x[[signal]][[sync]]))
    stream = x[[signal]][[sync]][[streamKey]]
  } else {
    if(!streamKey %in% c("s1","s2","time")) stop ("sync none requires streamkey == s1 or s2 or time")
    stream = x[[signal]][[streamKey]]
  }
  if( artefact.rm ){
    if(length(stream)!=length(x[[signal]]$valid)) stop("artefact.rm temporarily requires that stream has the same frequency of valid")
    stream[!x[[signal]]$valid]=NA
  }
  # seleziona il dataframe della categoria e crea un oggetto per ciascun livello
  cate = x[[category]]
  if(!is.factor(cate[[groupIndex]])) stop("groupIndex should be a factor column in the epoch table")
  resList = list()
  for(lev in levels(cate[[groupIndex]])){
    # dur = sum((cate[cate[[groupIndex]]==lev,"end"] - cate[cate[[groupIndex]]==lev,"start"] )*frequency(stream))
    if(mergeEpochs)
      resList[[lev]]=numeric()
    else
      resList[[lev]]=list()
  }
  names(resList) = levels(cate[[groupIndex]])

  remStream = stream
  for(i in 1:nrow(cate)){
    if(!is.na(cate[[groupIndex]][i]) && !is.null(cate[[groupIndex]][i])) {
      if(cate$start[i]>=end(stream)[1]){
        warning("In session ", dyadId(x),"-",sessionId(x), ", start of window ",i,": was equal to or beyond the stream end.", call.=F)
        lres[[i]] = NA
      } else { #if start is ok
        if(cate$end[i] > end(stream)[1] ){
          warning("In session ", dyadId(x),"-",sessionId(x), ", end of window ",i,": was reduced to the stream end.", call.=F)
          cate$end[i]= end(stream)[1]
        }
        # rimuovi la finestra dal segnale residuo:
        window(remStream, start = cate$start[i], end = cate$end[i]) <- NA
        
        #aggiungi la finestra al vettore di resList corrispondente al livello di groupIndex
        win = as.numeric(window(stream, start = cate$start[i], end = cate$end[i]))
        if(mergeEpochs)
          resList[[cate[[groupIndex]][i]]] = c(resList[[cate[[groupIndex]][i]]], win)
        else{
          resList[[cate[[groupIndex]][i]]] = c(resList[[cate[[groupIndex]][i]]], list(win))
          names(resList[[cate[[groupIndex]][i]]])[length(resList[[cate[[groupIndex]][i]]])] = paste0(dyadId(x),sessionId(x),"|",cate$start[i], "-",cate$end[i])
          }
      }
    }
  }
  resList[["remaining"]] = na.omit(as.numeric(remStream))
  
  #save object
  x[[signal]][[resName]] = resList 
  x
}




#### rinomina in "extractEpochs"
#### cambia epochStreamName in sync="PMBest", streamKey="sync", category="PACS/IM"
#' @export
#'
catExtractLong = function(experiment, signal="SC", epochStreamName="IM_PmdevSync", by, FUN = mean, ...){
  if(!missing("by")) stop("by is not implemented yet.")
  UseMethod("catExtractLong",experiment) 
}
#' @export
catExtractLong.DyadExperiment = function(experiment, signal, epochStreamName, by, FUN, ...){
  #check names
  keepNames = unique(unlist(lapply (experiment, function(session){
    names(session[[signal]][[epochStreamName]])
  })))
  keepNames = sort(keepNames)
  #controlla se mergeEpochs era T o F
  mergeEpochs = !any(unlist(lapply (experiment, function(session){
    lapply(session[[signal]][[epochStreamName]],is.list)
  })))
  #instanzia una lista vuota con keepNames elementi
  res = vector("list", length(keepNames))
  names(res) = keepNames

  for(session in experiment){
    if(!is.null(session[[signal]][[epochStreamName]])){
      sesNames= names(session[[signal]][[epochStreamName]])
      for(sesName in sesNames){
        newElement = session[[signal]][[epochStreamName]][[sesName]]
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

# questa funzione al momento è ibrida e incompleta.
# una prima parte serve ad estrarre delle finestre random dal remaining
# allo scopo di salvarle nel EXPERIMENT
# la seconda invece ripete la procedura 1000 volte allo scopo di fare un
# test di permutazione

# probabilmente ha senso tenere solo la seconda dal momento che la prima
# porta ad un database gigante. Eventualmente la prima si può astrarre con funzioni
# diverse dalla media.

#' Title 
#' Use either mimic or mean.duration and sd.duration
#'
#' @param mimic 
#' @param n if missing and mimic is used, the same number of the mimicked category is used
#' @param mean.duration 
#' @param sd.duration 
#' @param from 
#' @param experiment 
#' @param signal 
#' @param category 
#' @param sync 
#' @param streamKey 
# 
# randomEpochs = function(mimic, n, mean.duration, sd.duration, from="remaining", experiment, signal= "SC", category,  sync="PMBest", streamKey="sync"){
#   ##debug
#   mimic = "3"
#   n
#   mean.duration
#   sd.duration
#   from="remaining"
#   experiment = d3
#   signal= "SC"
#   sync="PMdev"
#   streamKey="sync"
#   category = "IM"
#   from="remaining"
#   #-----------------------
#   mimicMode = match.arg(mimicMode)
#   resName = paste0(c(toupper(category), "_", totitle(c(sync,streamKey))),collapse = "")
#   ex = catExtractLong(experiment, signal=signal, epochStream=resName)
#   ex2 = list(IM=do.call(c,ex[1:3]))
#   ex2$remaining = do.call(c,ex$remaining)
#   dur = lapply(ex, function(x)sapply(x,length))
#   
#   if(!from %in% names(ex)) stop ("from: ",from,"was not found in experiment")
#   
# 
#   if(!missing(mimic)) {
#     if(!missing(mean.duration) || !missing(sd.duration) ) stop ("specify either mimic or duration")
#     if(!mimic %in% names(ex)) stop ("mimic: ",mimic,"was not found in experiment")
#     if(missing(n)){
#       n = length(dur[[mimic]])
#       durList = dur[[mimic]]
#     } else { #se mimic, e 'n' è missing entra in strict mode
#       # e copia esattamente la durata delle finestre originali
#       mean.duration = mean(dur[[mimic]],na.rm=T)
#       sd.duration = sd(dur[[mimic]],na.rm=T)
#     }
#   }
#   if(!exists("durList")){ #a patto che non siamo in strict-mode
#     durList = c()
#     while(length(durList)<n) {
#       durList = c(durList, rnorm(n, mean.duration, sd.duration))
#       durList = durList[durList>0]
#       durList = sample(durList, min(length(durList),n), replace = F)
#     }
#     durList = round(durList)
#   }
#   ##ora seleziona le finestre
#   lfrom = do.call(c,ex[[from]]) #collassa tutti i remaining (o altro from) in un unica serie
#   newEpochs = list()
#   for(i in 1:n){
#     start = trunc(runif(1,1,length(lfrom)-max(durList)-1))
#     end = start + durList[i]
#     newEpochs[[i]] = lfrom[start:end]
#   }
#   boxplot(sapply(newEpochs,median),sapply(ex$`3`,median),ylim=c(0,1))
#   
#   #e se invece confrontassi la media reale con 1000 di quelle random?
#   n = length(ex2$IM)
#   dur2 = lapply(ex2, function(x)sapply(x,length))
#   durList = dur2$IM
#   lfrom = ex2$remaining
#   newEpochsMean = numeric(1000)
#   for(k in 1:1000){
#     newEpochs = list()
#     for(i in 1:n){
#       start = trunc(runif(1,1,length(lfrom)-max(durList)-1))
#       end = start + durList[i]
#       newEpochs[[i]] = lfrom[start:end]
#     }
#     newEpochsMean[k] = mean(sapply(newEpochs,median))
#   }
#   boxplot(newEpochsMean,sapply(ex2$IM,median),ylim=c(0,1))
#   length(newEpochsMean[newEpochsMean<mean(sapply(ex$`3`,median),na.rm=T)])
#   
#   
#   
# }
