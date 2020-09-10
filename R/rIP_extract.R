###################################

## NB
## cose con "apply" resitutiscono l'oggetto stesso trasformato (forse cambia nome? non è coerente con apply, lapply?)
## cose con "extract" restituiscono un oggetto diverso

#' Title
#' epochStreamApply cycles through all sessions in an experiment. For each session extracts a stream, then cuts it
#' according to categories epochs, and applies a function to each epoch.
#' 
#' @param x 
#' @param FUN 
#' @param signal 
#' @param sync either "none" or the name of a DyadSignal component, such as CCFBest or PMBest
#' @param streamKey 
#' @param category 
#' @param column 
#' @param ... arguments passed to FUN
#' 
#' @return the updated DyadExperiment, or single DyadSession object.
#' @export
#'
#' @examples
epochStreamApply = function(x, FUN, signal= "SC",sync = "none", streamKey=c("s1","s2","time","sync","lag","zero"),
                            category="PACS",... ){
  UseMethod("epochStreamApply", x)
}

#' @export
epochStreamApply.DyadExperiment = function(x, FUN, signal, sync, streamKey, category, artefact.rm = T, ...  ){
  FUNname = as.character(substitute(FUN))
  res = Map(function(session, nSession){
    prog(nSession, length(x))
    epochStreamApply(session, FUNname, signal, sync, streamKey, category, artefact.rm, ... )
  },x, seq_along(x))
  classAttr(res) = classAttr(x)
  res
}
#' @export
epochStreamApply.DyadSession = function(x, FUN, signal, sync , streamKey,
                                        category, artefact.rm = T, ...  ){
  if(!sync %in% names(x[[signal]])) stop ('sync should be one of:',names(x[[signal]]))
  
  
  if(is.character(FUN))
    FUNname = FUN
  else FUNname = as.character(substitute(FUN))
  FUN = match.fun(FUN)
  
  if(sync!="none"){
    if(!streamKey %in% names(x[[signal]][[sync]])) stop ('sync should be one of:',names(x[[signal]][[sync]]))
    stream = x[[signal]][[sync]][[streamKey]]
  } else {
    if(!streamKey %in% c("s1","s2","time")) stop ("sync none requires streamkey == s1 or s2 or time")
    stream = x[[signal]][[streamKey]]
    
  }
  if( artefact.rm ){
    stop("artefact.rm must be implemented with the new artefact data.frame architecture (in rIP_extract.R")
    if(length(stream)!=length(x[[signal]]$valid)) stop("artefact.rm temporarily requires that stream has the same frequency of valid")
    stream[!x[[signal]]$valid]=NA
  }
  
  
  cate = x[[category]]
  newName = paste(signal,sync,streamKey,FUNname,sep=".")
  # cate[,newName] = rep(NA,nrow(cate))
  lres = list()
  for(i in 1:nrow(cate)){
    if(cate$start[i]>=end(stream)[1]){
      warning("In session ", dyadId(x),"-",sessionId(x), ", start of window ",i,": was beyond the stream end.", call.=F)
      lres[[i]] = NA
    } else { #if start is ok
      if(cate$end[i] > end(stream)[1] ){
        message("In session ", dyadId(x),"-",sessionId(x), ", end of window ",i,": was reduced to the stream end.", call.=F)
        cate$end[i]= end(stream)[1]
      }
      win = window(stream, start = cate$start[i], end = cate$end[i])
      lres[[i]] = FUN(win, ...)
    }
    
  }
  funcols = names(lres[[1]])
  rexp2 = data.frame(do.call(rbind,lres))
  if(is.null(funcols)) colnames(rexp2) = FUNname
  else colnames(rexp2) = funcols
  colnames(rexp2) = paste(sync,streamKey,colnames(rexp2),sep = ".")
  
  cate = cbind(cate,rexp2)
  x[[category]]=cate
  x
}

#' Extracts summary data from categorial windows
#' the function first extracts all occurrences of a given category, combines them in a
#' single data.frame, then aggregates it according to specified rules.
#'
#' @param experiment 
#' @param category 
#' @param by  a character vector of one or more between "session","dyad","group" and the category's column names
#' @param FUN a function to compute the summary statistics which can be applied to all data subsets
#' @param ... arguments passed to FUN
#'
#' @return
#' @export
#'
#' @examples
catExtract = function(experiment, category="PACS", by, FUN = mean, ...){
  UseMethod("catExtract",experiment) 
}
#' @export
catExtract.DyadExperiment = function(experiment, category, by, FUN = mean, ...){
  #check names
  checkNames = unique(unlist(lapply (experiment, function(session){
    colnames(session[[category]])
  })))
  keepNames = checkNames
  for(iname in checkNames){
    if(!all(sapply(experiment,function(x){iname %in% colnames(x[[category]])})))
      keepNames = keepNames[keepNames!= iname]
  }
  if(!all(checkNames %in% keepNames)) message("'",paste(checkNames[!checkNames %in% keepNames], collapse="', '"),"' columns were not found in every session and were dropped")
  #rbind
  rexp = lapply(experiment, function(session){
    rses = session[[category]][,keepNames]
    rses$dyad = dyadId(session)
    rses$group = groupId(session)
    rses$session = sessionId(session)
    rses
  })
  res = do.call(rbind,rexp)
  
  #aggregate 
  if(missing(by)) by = NULL
  else if(!is.list(by)){
    by2 = match.arg(by, choices = c("session","dyad","group",keepNames), several.ok = T)
    if(!all(by %in% by2))message("'",paste(by[!by %in% by2], collapse="', '"),"' columns were not found and were ignored")
    by = lapply(by2, function(x){res[,x]})
    names(by) = by2
  }
  if(!is.null(by)){
    xNum = sapply(res, is.numeric)
    xNum[c("start","end","session","dyad","group")] = F
    res2 = aggregate(res[,xNum],by=by,FUN = FUN, ...)
    res3 = aggregate(res[,xNum],by=by,FUN = length)
    res2$n = res3[,length(res3)]
    res = res2
  }
  res
}


#' Extract summary data from a stream
#' The function iterates over each session of an experiment and applies a summarizing function on a given stream (e.g. mean).
#' Multiple output are possible (eg. using the quantile function)
#'
#' @param experiment 
#' @param sync 
#' @param streamKey 
#' @param FUN 
#' @param ... 
#' @param signal 
#'
#' @return
#' @export
#'
#' @examples
streamExtract<- function(experiment, signal, sync = "none", streamKey=c("s1","s2","time","sync","lag","zero"), FUN=mean, ...){
  UseMethod("streamExtract", experiment)
}
#' @export
streamExtract.DyadExperiment <- function(experiment, signal, sync, streamKey, FUN, ...){
  if(is.character(FUN))
    FUNname = FUN
  else FUNname = as.character(substitute(FUN))
  FUN = match.fun(FUN)
  
  rexp = Map(function(session,nsession){
    prog(nsession,length(experiment))
    if(sync!="none"){
      if(!streamKey %in% names(session[[signal]][[sync]])) stop ('sync should be one of:',names(session[[signal]][[sync]]))
      stream = session[[signal]][[sync]][[streamKey]]
    } else {
      if(!streamKey %in% c("s1","s2","time")) stop ("sync none requires streamkey == s1 or s2 or time")
      stream = session[[signal]][[streamKey]]
    }
    FUN(stream, ...)
    
  },experiment,seq_along(experiment))
  rexp2 = do.call(rbind,rexp)
  if(is.null(colnames(rexp2)))colnames(rexp2) = FUNname
  colnames(rexp2) = paste(sync,streamKey,colnames(rexp2),sep = ".")
  
  
  res = data.frame(sapply(experiment,groupId),
                   sapply(experiment,dyadId),
                   sapply(experiment,sessionId))
  names(res) = c("group","dyad","session")
  
  cbind(res, rexp2)
  
}
