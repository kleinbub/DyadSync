genericIO <- function (path,namefilt,idOrder,idSep, pairBind=F, ...){
  # if(missing(pairBind)) pairBind = F
  l = list(...)
  if(!file.exists(path)){stop("Selected file or directory does not exists")}
  patt = "\\.(txt|csv)$"
  if (!is.na(namefilt)){
    patt = paste0(namefilt,".*\\.(txt|csv)$")
  }
  nFiles=length(list.files(path, pattern=patt))
  if(nFiles==0){ #no files in the directory. Maybe a single file was provided?
    if(utils::file_test("-f",path)){ #check if file exists and not a directory
      filenames = path
      shortNames = basename(path)
      path=dirname(path)
    } else{
      stop("No [",patt,"] files were found in ",path,"\r\n\r\nFiles in path:\r\n",paste(list.files(path),"\r\n" ), call.=F)
      #stop("The input must be a file or a directory!\r\nTry to use choose.dir() as the first parameter, e.g.: meaMaster(choose.dir(), ....")
    }
  } else {
    filenames = list.files(path, pattern=patt, full.names=T)
    shortNames = list.files(path, pattern=patt, full.names=F)
  }
  nFiles = length(filenames)
  
  sapply(shortNames, cat, "\r\n")
  ####now tries to identify cases!
  
  #first check if separator exists
  if(any(!is.na(idOrder)))
    lapply(shortNames, function(name) {
      if(!grepl(idSep,tools::file_path_sans_ext(name))) stop("idSep: '",idSep, "' could not be found in file: ",name, call.=F)
    })
  
  #then build lists of groups, sessions, and id for each file
  #session
  idOrder = tolower(idOrder) #lowercase
  if(any(substr(idOrder,1,1)=="i",na.rm=T)){
    dyadIds = lapply(seq_along(shortNames), function(i) {
      unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="i")]
    })
  } else dyadIds = as.list(lead0(seq_len(nFiles),width=ifelse(nchar(nFiles)<2 ,2,nchar(nFiles)) ))
  #group
  if(any(substr(idOrder,1,1)=="g",na.rm=T)){
    group = lapply(seq_along(shortNames), function(i) {
      unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="g")]
    })
  } else group = as.list(rep("all",nFiles))
  #role
  if(any(substr(idOrder,1,1)=="r",na.rm=T)){
    if(!pairBind) stop("Role identifier was specified, but pairBind is false.")
    role = lapply(seq_along(shortNames), function(i) {
      unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="r")]
    })
    if(length(unique(role)) > 2 ) stop(paste(""))
  } else if(pairBind) {
    stop("pairBind = TRUE but no role identifier selected")
  } else role = as.list(rep("dyad",nFiles))
  #session
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
    
    nCheck = do.call(paste0, list(group,dyadIds, sess,role) ) #paste together group,id and session to obtain unique session identifier
    
    if(length(unique(nCheck)) != length(nCheck))
      warning("Two equal session identifier were found, please check that the session order/names are correct:\r\n",paste0("\r\n",shortNames,"\t",nCheck),call. = F)
    #now check if the session are in progressive order
    deltaSess = mapply(function(x,i){
      if(pairBind){
        if(sum(duplicated(x)) != length(unique(x))) stop ("Uncomplete dyad ", i,":\r\n", paste(x," ") )
        x = x[seq(1,length(x)-1,by=2)]
      }
      ifelse(length(x)>1,diff(x),1) #if there are multiple sessions with the same id, check that they are in progressive order
    }, split(unlist(sess),unlist(dyadIds)), unique(unlist(dyadIds) ) )
    if(any(deltaSess<=0))
      warning("The sessions may not be in sequential order, please check file names:\r\n",paste0("\r\n",shortNames,"\t",nCheck),call. = F)
  }
  
  
  ndyads = ifelse(pairBind, length(filenames)/2,length(filenames))
  nFiles = length(filenames)
  cat("\r\nReading",ndyads,"dyads\r\n")
  options = list(...)
  #ugly stuff to set read.csv options
  if("skip" %in% names(options)){
    if(length(options$skip) ==1) options$skip=rep(options[["skip"]],nFiles)
    if(length(options$skip) !=nFiles) stop("skip must be provided for each file")
    skipRow = options$skip
    options$skip = NULL
  } else skipRow =rep(0,nFiles)
  lf <- mapply(function(x,iFile) {  prog(iFile,nFiles); do.call(read.table,c(list(x, skip=skipRow[iFile]), options)) },filenames,seq_along(filenames),SIMPLIFY = F )
  if(ncol(lf[[1]])==1) {print(str(lf[[1]]));stop("Import failed. Check sep?")}
  return(list("lf"=lf,
              "sess"=sess,
              "dyadIds" = dyadIds,
              "group" = group,
              "role" = role,
              "filenames" = filenames,
              "shortNames" = shortNames,
              "ndyads" = ndyads,
              "nFiles" = nFiles
              ))
}