##readCategories v1.0 ############################################################
#this *UNTESTED* function should read any kind of savage csv format, and build
#a DyadCategories experiment.
#removeSec is a vector that specifies how many seconds to subtract from
#start and end columns of each file, useful to align categories to other signals.
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
readCategories = function(dirPath, startCol, endCol, catName, namefilt = NA,  removeSec = NULL, 
                          idOrder= c("dyadId","session","group"), #the order of identifiers in the filenames
                          idSep = "_", #the separator of identifiers in the filename
                          ... #additional options to be passed to read.table
){
  ###debug#####
  # dirPath = "A:\\OneDrive - Università degli Studi di Padova\\__Synchro Analysis\\__RLab\\longit_final\\PACS Talia\\C.C._seduta3_categorie-minutaggi_.csv"
  # dirPath = paste0(data_d,"\\giugno_2018")
  # options = list(sep=";")
  # startCol=1
  # endCol=2
  # removeSec = 277
  # sep=";"
  ##############
  
  
  iStep=0
  if(!file.exists(dirPath)){stop("Selected file or directory does not exists")}
  patt = "\\.(txt|csv)$"
  if (!is.na(namefilt)){
    patt = paste0(namefilt,".*\\.(txt|csv)$")
  }
  nFiles=length(list.files(dirPath, pattern=patt))
  if(nFiles==0){ #no files in the directory. Maybe a single file was provided?
    if(utils::file_test("-f",dirPath)){ #check if file exists and not a directory
      filenames = dirPath
      shortNames = basename(dirPath)
      dirPath=dirname(dirPath)
    } else{
      stop("No [",patt,"] files were found in ",dirPath,"\r\n\r\nFiles in path:\r\n",paste(list.files(dirPath),"\r\n" ), call.=F)
      #stop("The input must be a file or a directory!\r\nTry to use choose.dir() as the first parameter, e.g.: meaMaster(choose.dir(), ....")
    }
  } else {
    filenames = list.files(dirPath, pattern=patt, full.names=T)
    shortNames = list.files(dirPath, pattern=patt, full.names=F)
  }
  nFiles = length(filenames)
  
  
  ####now tries to identify cases!
  
  #first check if separator exists
  if(any(!is.na(idOrder)))
    lapply(shortNames, function(name) {
      if(!grepl(idSep,tools::file_path_sans_ext(name))) stop("idSep: '",idSep, "' could not be found in file: ",name, call.=F)
    })
  
  #then build lists of groups, sessions, and id for each file
  idOrder = tolower(idOrder) #lowercase
  if(any(substr(idOrder,1,1)=="i",na.rm=T)){
    dyadIds = lapply(seq_along(shortNames), function(i) {
      unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="i")]
    })
  } else dyadIds = as.list(lead0(seq_len(nFiles),width=ifelse(nchar(nFiles)<2 ,2,nchar(nFiles)) ))
  if(any(substr(idOrder,1,1)=="g",na.rm=T)){
    group = lapply(seq_along(shortNames), function(i) {
      unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="g")]
    })
  } else group = as.list(rep("all",nFiles))
  if(any(substr(idOrder,1,1)=="s",na.rm=T)){
    sess = lapply(seq_along(shortNames), function(i) {
      ax = unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="s")]
      x = as.numeric(gsub("[[:alpha:]]","", ax))
      if(is.na(x)){
        warning("No numeric information was found for session identifier '", ax, "' in signal ",shortNames[i],". Please check filenames and idOrder and idSep arguments:\r\n", call.=F)
        ax
      } else x
    })
  } else sess = as.list(rep("01",nFiles))
  
  
  #if sessions are specified, check their order
  if(any(substr(idOrder,1,1)=="s",na.rm=T) ){
    #now to check if the filenames have repetitions
    
    nCheck = do.call(paste0, list(group,dyadIds, sess) ) #paste together group,id and session to obtain unique session identifier
    if(length(unique(nCheck)) != length(nCheck))
      warning("Two equal session identifier were found, please check that the session order/names are correct:\r\n",paste0("\r\n",shortNames,"\t",nCheck),call. = F)
    
    #now check if the session are in progressive order
    deltaSess = sapply(split(unlist(sess),unlist(dyadIds)), function(x) ifelse(length(x)>1,diff(x),1)) #if there are multiple sessions with the same id, check that they are in progressive order
    if(any(deltaSess<=0))
      warning("The sessions may not be in sequential order, please check file names:\r\n",paste0("\r\n",shortNames,"\t",nCheck),call. = F)
  }
  
  
  iStep = iStep +1
  nFiles = length(filenames)
  cat("\r\nSTEP",iStep,"| Reading",nFiles,"dyads\r\n")
  options = list(...)
  if("skip" %in% names(options)){
    if(length(options$skip) ==1) options$skip=rep(options[["skip"]],nFiles)
    if(length(options$skip) !=nFiles) stop("skip must be provided for each file")
    skipRow = options$skip
    options$skip = NULL
  } else skipRow =rep(0,nFiles)
  lf <- mapply(function(x,iFile) {  prog(iFile,nFiles); do.call(read.csv,c(list(x, skip=skipRow[iFile]), stringsAsFactors = F, colClasses = "character",options)) },filenames,seq_along(filenames),SIMPLIFY = F )
  if(ncol(lf[[1]])==1) {print(str(lf[[1]]));stop("Import failed. Check sep?")}
  
  
  
  if(!is.null(removeSec) && length(removeSec)!=length(lf)) stop("removeSec must be defined for every file (n=",length(removeSec),")")
  
  cat("File name","\t","seconds removed\r\n")
  listCat = Map(function(file,iFile){
    cat(shortNames[iFile],"\t",removeSec[iFile],"\r\n")
    
    #remove na shit
    file[file==""] = NA
    file = file[rowSums(is.na(file)) != ncol(file),]
    file = file[,colSums(is.na(file)) != nrow(file)]
    file[is.na(file)] = "NA" #use character NA to keep smooth subsequent analyses
    #convert shitty time formats to seconds
    file[[startCol]] = timeMaster(file[[startCol]],-removeSec[iFile],out = "sec")
    file[[endCol]]   = timeMaster(file[[endCol]],-removeSec[iFile],out="sec")
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
    
    #convert remaining to factors
    data.frame("start"=file[[startCol]], "end" =file[[endCol]], "delta"=deltaSec, as.list(file[-c(startCol,endCol)]))
  },lf,seq_along(lf))
  
  
  experiment = DyadExperiment(dirPath,Map(function(session,nSession){
    #for each type of signal, add a new DyadSignal to the present DyadSession
    catList = list(
      DyadCategory(name=catName, data=session )
    )
    ses = DyadSession(groupId = group[[nSession]],sessionId= sess[[nSession]], dyadId = dyadIds[[nSession]],catList=catList,s1Name=NA,s2Name = NA, fileName = shortNames[[nSession]])
    names(ses$categ) = catName
    return(ses)
    
  },listCat,seq_along(listCat)))
  
  names(experiment) = paste(dyadIds,sess,sep="_")
  return(experiment)
}
