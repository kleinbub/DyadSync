

#' Title
#'
#' @param dirPath 
#' @param startCol 
#' @param endCol 
#' @param catName 
#' @param namefilt 
#' @param removeSec 
#' @param idOrder 
#' @param idSep 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
readCategories = function(path,
                          startCol, endCol, catName, namefilt = NA, 
                          removeSec = 0, 
                          idOrder= c("id","session","group"), #the order of identifiers in the filenames
                          idSep = "_", #the separator of identifiers in the filename
                          ... #additional options to be passed to read.table
){
  l = list(...) # l = list(sep=";")
  l = c(list(path,namefilt,idOrder,idSep),l)
  l$stringsAsFactors = F
  if(is.null(l$colClasses)) l$colClasses = "character"
  if(is.null(l$header)) l$header = TRUE
  
  imp = do.call(genericIO, l)
  lf= imp$lf
  sess= imp$sess
  dyadIds = imp$dyadIds
  group = imp$group
  role = imp$role
  filenames = imp$filenames
  shortNames = imp$shortNames
  ndyads = imp$ndyads
  nFiles = imp$nFiles
  
  if(length(removeSec)==1) {message("removeSec = ",removeSec," was used for all files")
    removeSec = rep(removeSec, length(lf))
  }  else if(length(removeSec)!=length(lf)) stop("removeSec must be defined for every file (n=",length(lf),")")

  cat("File name","\t","seconds removed\r\n")
  listCat = Map(function(file,iFile){
    cat(shortNames[iFile],"\t",removeSec[iFile],"\r\n")
    
    #remove na shit
    file[file==""] = NA
    file = file[rowSums(is.na(file)) != ncol(file),]
    file = file[,colSums(is.na(file)) != nrow(file)]
    file[is.na(file)] = "NA" #use character NA to keep smooth subsequent analyses
    #convert shitty time formats to seconds
    file[[startCol]] = timeMaster(file[[startCol]],out = "sec",add = -removeSec[iFile])
    file[[endCol]]   = timeMaster(file[[endCol]],out="sec",add = -removeSec[iFile])
    if(any(is.na(file[[startCol]]))) stop("NA was found in startCol on line ", which(is.na(file[[startCol]])))
    if(any(is.na(file[[endCol]])))   stop("NA was found in startCol on line ", which(is.na(file[[endCol]])))
    
    deltaSec = file[[endCol]] - file[[startCol]]
    deltaSec[deltaSec==0] =1
    if(any(file[[endCol]]<0))   stop("After removeSec negative times were found in 'end' column in file "  ,shortNames[iFile])
    if(any(file[[startCol]]<0)) stop("After removeSec negative times were found in 'start' column in file ",shortNames[iFile])
    
    if(length(deltaSec[deltaSec<0])>0) warning("negative durations spotted in file ",shortNames[iFile])
    
    
    
    #trim leading and ending whitespaces
    for(i in 1:ncol(file)){
      if(is.character(file[[i]])){
        file[[i]] = gsub("^\\s+|\\s+$", "",  file[[i]])
      }
    }
    
    
    res = data.frame("start"=file[[startCol]], "end" =file[[endCol]], "delta"=deltaSec)
    #all other columns should be characters if colClasses was not overridden
    res = cbind(res, file[-c(startCol,endCol)])
    #check it and finally apply type.convert on all other columns
    k = list(...)
    if(is.null(k$colClasses)){ 
      res[which(sapply(res,is.character))] = lapply(res[which(sapply(res,is.character))],type.convert)
    }
    res
  },lf,seq_along(lf))
  
  
  experiment = DyadExperiment(path,Map(function(session,nSession){
    #for each type of signal, add a new DyadSignal to the present DyadSession
    signalList = list(
      DyadCategory(name=catName, data=session )
    )
    ses = DyadSession(groupId = group[[nSession]], sessionId= sess[[nSession]], dyadId = dyadIds[[nSession]],
                      signalList=signalList, s1Name = NA, s2Name = NA, fileName = shortNames[[nSession]] )
    names(ses) = catName
    return(ses)
    
  },listCat,seq_along(listCat)))
  
  names(experiment) = paste(dyadIds,sess,sep="_")
  return(experiment)
}
