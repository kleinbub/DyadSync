##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  
## R class dyadSignal
# Main class definitions 
############################################################################################################ 
## Credits
# Author: Johann R. Kleinbub
# Contact: johann.kleinbub@gmail.com
############################################################################################################ 
##Development notes regarding streams
# L'idea di base è che ogni segnale è uno stream. Questo già avviene per i dati grezzi, ma non ancora per i ccf
# o per i dati categoriali. la matrice con tutti i lag, andrebbe abbandonata, o quanto meno salvata a parte e 
# bestCCF e bestLAG salvati come stream. Gli stream infatti sono time-series con un sacco di metadati utili
# che rendono molto più semplice avere sempre tutte le info a portata di mano

## STRUTTURA v3
## DyadExperiment               | <- Experiment è semplicemente il contenitore di tutti i dati
##    $DyadSession              | <- la session è l'unità logica, contiene nomi partecipanti e id
##        $DyadSignal "SC"      | <- signal è il contenitore di tutti i dati e analisi di un tipo di segnale fisio
##            $DyadStream s1    | <- i diversi dati sono in forma di stream, ossia serie temporali con 
##            $DyadStream s2    |    più metadati, tra cui elementi grafici che forse sono da eliminare
##            $artefact         | <- un data.frame contenente le epoche (start end in secondi) da escludere
##            $CCFBest          | <- contenitore di analisi sincro basato sulle windowed x-cors
##              $rats sync      | <- il vecchio BestCCF
##              $rats lag       | <- il vecchio BestLag
##              $df table       | <- la vecchia ccfMatrix
##            $AMICo            | <- contenitore di analisi sincro basato sul Peak Matching
##              $rats sync      |
##              $rats lag       |
##              $df xBest       |
##        $DyadSignal "PPG"     |
##            $rats s1          |
##            $rats s2          |
##            $rats ...         |
##        $DyadCategory "PACS"  | <- oltre ai segnali, una sessione contiene le categorie, che sono dataframe contenenti
##                              |    finestre temporali di interesse


### cose che deve contenere CCF:
##BestCCF PROCEDURE
# matrice di xCor di dimensioni nLag x nWin
# stream di bestCCF a BestLag
# stream di bestLag
# valori riassuntivi bestCCF
## PPSync procedure
# tabella xbest 
# stream di sync secondo ppBest
# stream di lag secondo ppBest
# valori riassuntivi ppSync
# valori di zSync per differenti metodi

## proposta 1
#  oggetto di classe DyadSync come contenitore generico di analisi di sincronizzazinoe
## -contiene "nome" o "tipo" es: ccf, best procedure, PP procedure,
## -ha dei metodi generalizzati per questa classe es toStream() che indipendentemente dal tipo
##  restituisce uno stream interpolato alla frequenza richiesta.
## -però è brutto generalizzare i metodi sull'attributo "nome"
## proposta 2
## - creare classi diverse e definire i metodi per le diverse classi
## -ccfMatrix, vectorBestLag, peakpicking best lag,

## ----------------------------------------------------------




### DYADEXPERIMENT ###########################################
## list of DyadSessions with attributes:
##   -name
##   -class: list DyadExperiment
DyadExperiment = function(name, sessionList){
  exp =  sessionList
  attr(exp, "name") = name
  class(exp) = "DyadExperiment"
  return(exp)
}
#' @export
is.DyadExperiment = function(x) inherits(x,"DyadExperiment") && length(x)

#' Extract or replace parts of a dyadexperiment
#'
#' @param x 
#' @param i 
#' @param name The name of the new experiment
#'
#' @return
#' @export
#'
#' @examples
"[.DyadExperiment" = function(x,i,name=NA){
  y = .subset(x, i)
  if(is.na(name)) name = paste0(paste0(name(x),collapse="-" ),"-redux")
  DyadExperiment(name,y)
}

# c.DyadExperiment sostituisce experimentMerge() e permette di unire
# i segnali di esperimenti con la stessa struttura di sessioni.
#' @export
c.DyadExperiment = function (...){
  l = list(...)
  
  #######################
  # l = list(mea,sc)
  if(length(l)==1) return(l)#return(l[[1]])
  
  
  
  #c deve unire tra di loro i segnali che hanno lo stesso sessionId, Group, dyadID
  #invece deve incollare uno dopo l'altro le diadi che non overlappano
  #se per lo stesso ID ci sono due volte lo stesso segnale, stop.
  # comp = lapply(l,sapply,function(session){paste(groupId(session),dyadId(session),sessionId(session),sep="_")})
  # comp = lapply(comp,tolower)
  # 
  if(length(unique(sapply(l,class)))>1) stop("Only 'DyadExperiment' objects can be combined together")
  
  newEX = unclass(l[[1]])
  newEXsessions = attr(newEX,"names")
  fancynames = paste(newEXsessions, sapply(lapply(newEX, names), paste,collapse="_"),sep="_")
  report = data.frame("FINAL"=fancynames, "exp1"= fancynames)
  joined = added = 0
  
  for(y in (2):length(l)){ #per tutti gli esperimenti successivi
    #to compare
    ADD = l[[y]]
    ADDsessions = attr(ADD,"names")
    
    report = cbind(report,data.frame(temp=NA))
    colnames(report)[y+1] = paste0("exp",y)
    
    
    for(s in 1:length(ADD)){#per ciascuna sessione di ADD
      # print(s)
      toJoin = which(newEXsessions ==ADDsessions[s])
      if(length(toJoin) > 1) stop("experiment ", y, " had more than one session of ",ADDsessions[s])
      if(length(toJoin) ==1){ #if there was a match
        # check that signals are not 
        if(any(names(ADD[[ADDsessions[s]]]) %in% names(newEX[[newEXsessions[toJoin]]]) )) stop ("session ",ADDsessions[s]," in experiment 1 and ",y," had the same signal")
        else {
          report[toJoin,1] = paste(report[toJoin,1], paste(names(ADD[[s]]),collapse = "_"),sep="_") 
          report[toJoin,y+1] = paste(ADDsessions[s], paste(names(ADD[[s]]),collapse = "_"),sep="_") 
          joined = joined +1
          newEX[[toJoin]] = c(newEX[[toJoin]],ADD[[s]]) #aggiungi all'exp originale
        }
      } else { #the session could not be joined, so add it! 
        newEX = c(newEX,ADD[s]) 
        report[nrow(report)+1,1] = paste(ADDsessions[s], paste(names(ADD[[s]]),collapse = "_"),sep="_")
        report[nrow(report),y+1] = paste(ADDsessions[s], paste(names(ADD[[s]]),collapse = "_"),sep="_")
        added = added+1
      }
      #refresh names of newEX
      newEXsessions = attr(newEX,"names")
      
    }
  }
  
  print(report)
  cat ("\r\nMerge successful.",added+length(l[[1]]),"sessions were added, and ",joined," were joined to existing sessions. The final DyadExperiment consists of", length(newEX),"unique sessions.")
  
  DyadExperiment(sapply(l,name),newEX)
}

### DYADSESSION ##############################################
## list of signals and categ with attributes:
##   -name
##   -sessionId
##   -dyadId
##   -groupId
##   -s1Name
##   -s2Name
##   -fileName
##   -class: list DyadSession
DyadSession = function(groupId,sessionId,dyadId, signalList=NULL, s1Name,s2Name,fileName){
  if(length(unique(sapply(signalList, "s1Name")))>1)stop("multiple s1Names are not supported")
  if(length(unique(sapply(signalList, "s2Name")))>1)stop("multiple s2Names are not supported")
  
  x = signalList
  
  if(!is.null(signalList)) names(x) = lapply(x, name)
  attributes(x) = c(attributes(x),list(
    name = paste(groupId,dyadId,lead0(sessionId),sep="_"),
    sessionId = sessionId,
    dyadId = dyadId,
    groupId = groupId,
    s1Name = s1Name,
    s2Name = s2Name,
    fileName = fileName
  ))
  class(x) = "DyadSession"
  return(x)
}
#' @export
is.DyadSession = function(x) inherits(x,"DyadSession") && length(x)

#' @export
c.DyadSession = function(...){
  l = list(...);
  x = l[[1]]
  fileNames = sapply(l, attr, "fileName")
  #check on different filenames (bad)
  # if(length(unique(sapply(l,name))) >1 ) stop("Can't combine sessions with different names:\r\n",paste(unlist(sapply(l, attr, "fileName")),collapse="\r\n"), call.=F)
  #check on same signal names (bad)
  if(length(unique(sapply(l,names)))<length(sapply(l,names))) stop("Can't combine sessions containing the same signal. Use selectSignals() to extract only different signals before merging", call.=F)
  #check on different s1 s2 names (bad)
  if(any(na.omit(sapply(l,s1Name)) %in% na.omit(sapply(l,s2Name)))) warning("Sessions contain different s1 and s2 names", call.=F)
  
  structure(NextMethod("c"),
            "name" = name(x),
            "sessionId" = sessionId(x),
            "dyadId" = dyadId(x),
            "groupId" = groupId(x),
            "s1Name" = na.omit(sapply(l,s1Name)), #even if one of the signals ha NA sNames
            "s2Name" = na.omit(sapply(l,s2Name)), #all values are selected
            "fileName" = fileNames,
            "class" ="DyadSession")
}


### DYADCATEGORY ###########################################
## data.frame with attributes:
##   -name
##   -class: list DyadCategory
DyadCategory = function(name, data){
  categ = data
  attributes(categ) = c(attributes(categ),list(name = name))
  class(categ) = append(class(categ),"DyadCategory")
  return(categ)
}
#' @export
is.DyadCategory = function(x) inherits(x,"DyadCategory") && length(x)

### DYADSIGNAL ###########################################
## list of DyadStreams, time, valid with attributes:
##   -SR
##   -filter 
##   -ccf 
##   -s1Name
##   -s2Name
##   -start
##   -end
##   -duration
##   -class: list DyadSession
DyadSignal = function(name="some signal",s1=NULL,s2=NULL,SR=NULL,
                      s1Name, s2Name,sessionId,dyadId,groupId){
  if(!is.rats(s1) || !is.rats(s2)) stop("both s1 and s2 must be rats")
  x = list(
    s1 = s1,
    s2 = s2,
    artefacts = data.frame("start"=c(),"end"=c())
  )
  attributes(x) = c(attributes(x),list(
    "name" = name,
    "sessionId" = sessionId,
    "dyadId" = dyadId,
    "groupId" = groupId,
    "s1Name" = s1Name,
    "s2Name" = s2Name,
    "SR" = SR,
    "filter" = "raw",
    "start" = start(s1), #start-end-duration of s1 and s2 are the same by design
    "end" = end(s1)
    
  ))
  class(x) = append(class(x),"DyadSignal")
  return(x)
} 
#' @export
is.DyadSignal = function(x) inherits(x,"DyadSignal") && length(x)

#' Time Windows
#'
#' @param x 
#' @param duration an alternative specification to end
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
window.DyadSignal = function(x, duration=NULL, ...){
  l = list(...)
  l = setArg("start",attr(x,"start"),l)
  if(!is.null(duration)){
    xend =  c(l[["start"]][1] + duration-1, frequency(x)) 
  } else xend = attr(x,"end")
  l = setArg("end",xend,l)
  #1. find all rats in x and nested objects
  #2. window them all
  #3. recreate meta-data
  my_rats = which(unlist(lapply(x,is.rats)))
  my_rats = c(my_rats,which(unlist(lapply(x,is.sync))))
  res = x
  for(i in my_rats){
    if(is.rats  (res[[i]])) res[[i]] = do.call("window",c(list(res[[i]]),l))
    if(is.sync(res[[i]])){
      my_rats2 = which(unlist(lapply(res[[i]],is.rats)))
      for(j in my_rats2){
        res[[i]][[j]] = do.call("window",c(list(res[[i]][[j]]),l))
      }
    }
  }
  
  classAttr(res) = classAttr(x)
  # attr(res,"duration") = length(res)/frequency(res)
  attr(res,"start") = start(res)
  attr(res,"end") = end(res)
  res
}


