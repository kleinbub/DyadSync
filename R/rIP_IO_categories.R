

#' Title
#'
#' @param startCol character. The column containing the start of each epoch.
#' In hh:mm:ss or mm:ss format, or just a number representing seconds (see timeMaster).
#' @param endCol character. The column containing the end of each epoch.
#' @param codingCol numeric vector. The column(s) containing categories. These
#' will be converted to factors, searching across all files for possible levels.
#' @param catName character. The name of the category
#' @param namefilt character. A name or part of name, to filter among all files in path.
#' @param removeSec numeric or numeric vector. number of seconds to be subtracted from the 
#' @param idOrder 
#' @param idSep 
#' @param ... 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
readCategories = function(path,
                          startCol, endCol,codingCol, catName, namefilt = NA, 
                          removeSec = 0, 
                          idOrder= c("id","session","group"), #the order of identifiers in the filenames
                          idSep = "_", #the separator of identifiers in the filename
                          ... #additional options to be passed to read.table
){
  # debug
  # path = "A:/OneDrive - Università degli Studi di Padova/__Ricerche/2019_04_Brasini/PIRS_brasini"
  # namefilt = ""
  # catName = "PIRS"
  # removeSec = 0
  # startCol =1
  # endCol = 2
  # sep = ";"
  # idOrder = c("id", "session")
  # idSep = "_"
  ##############
  
  l = list(...) # l = list(sep=";")
  l = c(list(path,namefilt,idOrder,idSep),l)
  l$stringsAsFactors = FALSE
  l$as.is = TRUE
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
    
    for(i in 1:length(file[[startCol]])) {
      result = tryCatch({
        a = timeMaster(file[[startCol]][i],out = "sec",add = -removeSec[iFile])
        b = timeMaster(file[[endCol]][i],out = "sec",add = -removeSec[iFile])
        if(b-a <0 )
          cat("--negative duration spotted at line:",i,"\r\n")
        if(a<0 )
          cat("--negative start time (after removeSec) spotted at line:",i,"\r\n")
        if(b<0 )
          cat("--negative end time (after removeSec) spotted at line:",i,"\r\n")
        if(i>1){ #check if start is after the previous end. Epochs should NOT overlap
          if(a<timeMaster(file[[endCol]][i-1],out = "sec",add = -removeSec[iFile]))
            cat("--overlapping epochs spotted at line:",i,"\r\n")
          }
        
      }, error = function(e) {
        cat("--invalid time value found at line:",i,"\r\n")
      }, finally = {
      })
    }
    
    
    
    file[[startCol]] = timeMaster(file[[startCol]],out = "sec",add = -removeSec[iFile])
    file[[endCol]]   = timeMaster(file[[endCol]],out="sec",add = -removeSec[iFile])

    if(any(is.na(file[[startCol]]))) stop("NA was found in startCol on line ", which(is.na(file[[startCol]])))
    if(any(is.na(file[[endCol]])))   stop("NA was found in startCol on line ", which(is.na(file[[endCol]])))
    
    deltaSec = file[[endCol]] - file[[startCol]]
    deltaSec[deltaSec==0] =1
    # if(any(file[[endCol]]<0))   stop("After removeSec negative times were found in 'end' column in file "  ,shortNames[iFile])
    # if(any(file[[startCol]]<0)) stop("After removeSec negative times were found in 'start' column in file ",shortNames[iFile])
    
    #if(length(deltaSec[deltaSec<0])>0) warning("negative durations spotted in file ",shortNames[iFile],"\r\n",which(deltaSec[deltaSec<0],arr.ind = T))
    
    
    
    #trim leading and ending whitespaces
    for(i in 1:ncol(file)){
      if(is.character(file[[i]])){
        file[[i]] = gsub("^\\s+|\\s+$", "",  file[[i]])
      }
    }
    
    
    res = data.frame("start"=file[[startCol]], "end" =file[[endCol]], "delta"=deltaSec)
    #add the coding columns, which will be converted to factor later
    res = cbind(res, file[, codingCol])
    #add all other columns that should be characters if colClasses was not overridden
    res = cbind(res, file[-c(startCol,endCol,codingCol)])
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
  
  #the codingCols have been placed together starting from column 4
  newCodingCol = 4:(3+length(codingCol))
  
  #now for every coding, run across all files to search all possible levels, then transform to factor.
  for(cate in 1:length(newCodingCol)){
    xxx = c()
    yyy = c()
    for(i in 1:length(experiment)){
      xxx=c(xxx,experiment[[i]][[catName]][,newCodingCol[cate]])
      yyy=c(yyy,colnames(experiment[[i]][[catName]][newCodingCol[cate]]))
    }
    codes = sort(unique(xxx))
    titles = sort(unique(yyy))
    
    #are there uncoherent colnames?
    if(length(titles)>1){
      cat0("\r\n\r\nIn column ",codingCol[cate], " various names where found and coherced to the most common:\r\n")
      cat(titles)
      titles = names(table(yyy))[which.max(table(yyy))]
      for(i in 1:length(experiment)){
        colnames(experiment[[i]][[catName]][newCodingCol[cate]]) = titles
      }
    }
    cat0("\r\n\r\nFor coding '",titles,"', the following levels where detected:\r\n")
    table(xxx)
    cat0("\r\n --------------------")
    for(i in 1:length(experiment)){
      experiment[[i]][[catName]][,newCodingCol[cate]] = factor(experiment[[i]][[catName]][,newCodingCol[cate]],levels = codes )
    }
    
  }
  return(experiment)
}
