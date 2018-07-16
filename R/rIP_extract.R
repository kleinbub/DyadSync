###################################

## NB
## cose con "apply" resitutiscono l'oggetto stesso trasformato (forse cambia nome? non è coerente con apply, lapply?)
## cose con "extract" restituiscono un oggetto diverso

#' Title
#'
#' @param x 
#' @param FUN 
#' @param signal 
#' @param sync 
#' @param streamKey 
#' @param category 
#' @param column 
#' @param ... arguments passed to FUN
#'
#' @return
#' @export
#'
#' @examples
catApply = function(x, FUN, signal= "SC",sync = c("none",SYNC_CLASSES), streamKey=c("s1","s2","time","sync","lag","zero"),
                    category="PACS",... ){
  UseMethod("catApply", x)
}

#' @export
catApply.DyadExperiment = function(x, FUN, signal, sync, streamKey, category, ...  ){
  FUNname = as.character(substitute(FUN))
  res = Map(function(session, nSession){
    prog(nSession, length(x))
    catApply(session, FUNname, signal, sync, streamKey, category, ... )
  },x, seq_along(x))
  classAttr(res) = classAttr(x)
  res
}
#' @export
catApply.DyadSession = function(x, FUN, signal, sync = c("none",SYNC_CLASSES), streamKey=c("s1","s2","time","sync","lag","zero"),
                                category,...  ){
  sync = match.arg(sync)
  streamKey = match.arg(streamKey)

  if(is.character(FUN))
    FUNname = FUN
  else FUNname = as.character(substitute(FUN))
  FUN = match.fun(FUN)
  
  if(sync!="none"){
    stream = x[[signal]][[sync]][[streamKey]]
  } else {
    if(!streamKey %in% c("s1","s2","time")) stop ("sync none requires streamkey == s1 or s2 or time")
    stream = x[[signal]][[streamKey]]
  }
  
  cate = x[[category]]
  newName = paste(signal,sync,streamKey,FUNname,sep=".")
  cate[,newName] = rep(NA,nrow(cate))
  for(i in 1:nrow(cate)){
    if(end(stream)[1]>cate$end[i]){
      win = window(stream, start = cate$start[i], end = cate$end[i])
      cate[i,newName] = FUN(win,na.rm=T)#FUN(win,...)
    } else {
      warning("In session ", sessionId(session), ", end of window ",i,": was beyond the stream end.")
    }
  }
  x[[category]]=cate
  x
}

#' Extract summary data from categorial windows
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
      keepNames = checkNames[checkNames!= iname]
  }
  if(!all(checkNames %in% keepNames)) warning("'",paste(checkNames[!checkNames %in% keepNames], collapse="', '"),"' columns were not found in every session and were dropped")
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
    if(!all(by %in% by2))warning("'",paste(by[!by %in% by2], collapse="', '"),"' columns were not found and were ignored")
    by = lapply(by2, function(x){res[,x]})
    names(by) = by2
  }
  if(!is.null(by)){
    xNum = sapply(res, is.numeric)
    xNum[c("start","end") ] = F
    res = aggregate(res[,xNum],by=by,FUN = FUN, ...)
  }
  res
}


#' Title
#'
#' @param experiment 
#' @param sync 
#' @param streamKey 
#' @param by 
#' @param FUN 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
streamExtract<- function(experiment, signal, sync = c("none",SYNC_CLASSES), streamKey=c("s1","s2","time","sync","lag","zero"), FUN=mean, ...){
  UseMethod("streamExtract", experiment)
}
#' @export
streamExtract.DyadExperiment <- function(experiment, signal, sync = c("none",SYNC_CLASSES), streamKey=c("s1","s2","time","sync","lag","zero"), FUN=mean, ...){
  sync = match.arg(sync)
  streamKey = match.arg(streamKey)

  if(is.character(FUN))
    FUNname = FUN
  else FUNname = as.character(substitute(FUN))
  FUN = match.fun(FUN)

  rexp = Map(function(session,nsession){
    prog(nsession,length(experiment))
    if(sync!="none"){
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
 