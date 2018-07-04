##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  
## R class dyadSignal
# Main class definitions and libraries importer
#
############################################################################################################ 
## Changelog
# v2.0 integrato in rIP package
# v1.6 usa la nuova versione di CCF, non backward compatibile
# v1.4 a lot of goodies :3
# v1.3 attributes instead than lists. Not backward compatible
# v1.1 added "DyadStream" concept
# v1.0 stable
#
############################################################################################################ 
## ToDo *asap*
#
# -everything should become streams. and one should add as many streams as needed, es: many ccf with various settings.
# -every class should use attributes instead of lists, where appropriate
## Categorial Love
# -DyadCategory makes sense as a dataframe. but it needs adequate methods, for instance:
# -DyadCategory to stream (category, column) : crea versione interpolata. utile principalmente solo per plottare
# -byCat(signal, category, stream=c(patient, therapist, bestCCF, bestLag,lag0), FUN) : applica una funzione es la media ai diversi livelli
#       del factor di uno stream. a seconda che input sia experiment, session o singal applica adeguatamente.
# ----> queste cose sono parzialmente applicate in DyadCategory_BETA
#
## ToDo *later*
#
# REMEMBER: ####
# -"end" and "start" are complex. And strictly interrelated with frequency.
# -Very important! Use start = 0, otherwise things get exotheric.
# -Eg: a <- ts(runif(147), start=0, frequency = f)
# -if f = 10: end(a) == c(14, 7) means "14.7" seconds [ != tsp(a)[2L] remember!]
# -if f =  5: end(a) == c(29, 2) means "29 times full 5 items, and only the first 2 elements of the 30th cycle" (29x5 + 2 = 147)
#   since each element is 1/5 of a cycle (0.2s) the total time is 29.4 seconds
# -The code forces start = c(0,1)
# start ed end ragionano come se fossero mesi. quindi start(a) == 0 1, significa il "primo mese dell'anno zero". L'idea è che la
# frequenza abbia un significato, es 1 mese, 1 osservazione, quindi end 29 2 significa la seconda osservazione del 29 ciclo.
# siccome il mio ciclo sono i secondi, usare come frequenza 10 o 100 rappresenta decimi e centesimi di secondo, che è sensato.
# usare frequenza 5 significa 
# 
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

## STRUTTURA:
## DyadExperiment               | <- Experiment è semplicemente il contenitore di tutti i dati
##    $DyadSession              | <- la session è l'unità logica, contiene nomi partecipanti e id
##        $DyadSignal "SC"      | <- signal è il contenitore di tutti i dati e analisi di un tipo di segnale fisio
##            $DyadStream s1    | <- i diversi dati sono in forma di stream, ossia serie temporali con 
##            $DyadStream s2    |    più metadati, tra cui elementi grafici che forse sono da eliminare
##            $DyadStream ccf1  |
##            $DyadStream ccf2  |
##            $DyadStream ...   |
##            $CcfMatrix          | <- un vecchio contenitore della matrice di correlazione e i metadati con i setting della ccf
##                              |    devo trovare il modo di renderlo più generale per accomodare diverse forme di ccf
##        $DyadSignal "PPG"     |
##            $DyadStream s1    |
##            $DyadStream s2    |
##            $DyadStream ...   |
##        $DyadCategory "PACS"  | <- oltre ai segnali, una sessione contiene le categorie, che sono dataframe contenenti
##                              |    finestre temporali di interesse



### DYADEXPERIMENT ###########################################
## list of DyadSessions with attributes:
##   -name
##   -class: list DyadExperiment
DyadExperiment = function(name, sessionList){
  exp =  sessionList
  attributes(exp) = c(attributes(exp), list(name   = name))
  class(exp) = append(class(exp),"DyadExperiment")
  return(exp)
}
is.DyadExperiment = function(x) inherits(x,"DyadExperiment") && length(x)


### DYADSESSION ##############################################
## list of signals and categ with attributes:
##   -sessionId
##   -dyadId
##   -uid
##   -s1Name
##   -s2Name
##   -class: list DyadSession
DyadSession = function(sessionId,dyadId, signalList=NULL,catList=NULL, s1Name,s2Name){
  if(length(unique(sapply(signalList, "s1Name")))>1)stop("multiple s1Names are not supported")
  if(length(unique(sapply(signalList, "s2Name")))>1)stop("multiple s2Names are not supported")
  
  x = list(
    signals = signalList,
    categ = catList
  )
  if(!is.null(signalList)) names(x$signals) = lapply(x$signals, function(x) x$name)
  attributes(x) = c(attributes(x),list(
    sessionId = sessionId,
    dyadId = dyadId,
    uid = paste(sessionId,dyadId,sep="_"),
    s1Name = s1Name,
    s2Name = s2Name
  ))
  class(x) = append(class(x),"DyadSession")
  return(x)
}
is.DyadSession = function(x) inherits(x,"DyadSession") && length(x)


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
is.DyadCategory = function(x) inherits(x,"DyadCategory") && length(x)

### DYADSIGNAL ###########################################
## list of DyadStreams, CcfMatrix, and sampRate with attributes:
##   -sampRate
##   -filter 
##   -ccf 
##   -s1Name
##   -s2Name
##   -start
##   -end
##   -duration
##   -class: list DyadSession
DyadSignal = function(name="some signal",s1=NULL,s2=NULL,sampRate=NULL, s1Name, s2Name){
  x = list(
    name = name,
    s1 = DyadStream(stream = s1, name=s1Name, sampRate = sampRate, col = "deeppink3",  lty=1, lwd=2),
    s2 = DyadStream(stream = s2, name=s2Name, sampRate = sampRate, col = "dodgerblue3", lty=1, lwd=2),
    valid = DyadStream(stream = ts(rep(TRUE, length(s1)),start=start(s1),end= end(s1), frequency = frequency(s1) ), name="valid", sampRate = sampRate, col = "dodgerblue3", lty=1, lwd=2),
    sampRate = sampRate,
    ccf = NULL
  )
  attributes(x) = c(attributes(x),list(
    sampRate = sampRate,
    filter = "raw",
    ccf = FALSE,
    s1Name = s1Name,
    s2Name = s2Name,
    start = start(s1), #start-end-duration of s1 and s2 are the same
    end = end(s1),
    duration = trunc(length(s1)/sampRate)
  ))
  class(x) = append(class(x),"DyadSignal")
  return(x)
} 

is.DyadSignal = function(x) inherits(x,"DyadSignal") && length(x)

### CcfMatrix ###########################################
## list of of ccfmat, sampRate, and settings with attributes:
##   -class: list CcfMatrix

CcfMatrix = function(ccfmat, sampRate,
                     settings = list("lagSec"=NULL,"incSec"=NULL,"winSec"=NULL,"accelSec"=NULL,"weight"=NULL,"interpolated"=NULL)){
  CCF =  list(
    ccfmat   = ccfmat,
    sampRate = sampRate,
    settings = list(
      lagSec   = settings$lagSec,
      incSec   = settings$incSec,
      winSec   = settings$winSec,
      accelSec = settings$accelSec,
      weight = settings$weight,
      interpolated = settings$interpolated
    )
  )
  class(CCF) = append(class(CCF),"CcfMatrix")
  return(CCF)
  

  }

### DYADSTREAM ###########################################
## time-serie ts() with additional attributes:
##   -name
##   -settings 
##   -col 
##   -lty
##   -lwd
##   -tsp (inherited by ts)
##   -class: list DyadSession
DyadStream = function(stream, name, sampRate=frequency(stream), start=0, settings = list(), col=1, lty=1, lwd=1){
  if(!is(stream,"ts")){
    stream = ts(stream, start = start, frequency = sampRate)
    warning(paste0("Stream '",name,"' was coerced to ts starting at 0, with sampRate of: ",sampRate,"Hz, and duration of: ",end(stream)[1]," seconds."), call.=F)
  }
  attributes(stream) = c(attributes(stream),
                         list(
                            name = name,
                            duration = trunc(length(stream)/frequency(stream)),
                            settings = settings,
                            col=col,
                            lty=lty,
                            lwd = lwd
                        ))
  class(stream) = append("DyadStream",class(stream))
  return(stream)
}

is.DyadStream = function(x) inherits(x,"DyadStream") && length(x)

print.DyadStream = function (x, ...) {
  #x <- as.ts(x)
  Tsp <- tsp(x)
  if (is.null(Tsp)) {
    warning("series is corrupt, with no 'tsp' attribute")
    print(unclass(x), ...)
    return(invisible(x))
  }
  nn <- 1 + round((Tsp[2L] - Tsp[1L]) * Tsp[3L])
  if (NROW(x) != nn) {
    warning(gettextf("series is corrupt: length %d with 'tsp' implying %d", 
                     NROW(x), nn), domain = NA, call. = FALSE)
    calendar <- FALSE
  }
  fr.x <- frequency(x)
  if (fr.x != 1) 
    cat0("DyadStream '",attr(x,"name"),"':\nStart: ", start(x)[1]+((start(x)[2]-1)/frequency(x)), " seconds =",deparse(start(x)),
        "\nEnd =", end(x)[1]+(end(x)[2]/frequency(x)), " seconds =",deparse(end(x)) ,
        "\nFrequency = ", deparse(fr.x), 
        "\n")
  else cat0("DyadStream '",attr(x,"name"),"':\nStart: ", format(tsp(x)[1L]), 
           "\nEnd: ", format(tsp(x)[2L]), "\nFrequency: ", deparse(fr.x), 
           "\n")
  cat0("Duration: ",length(x), " samples\n")
  cat0("start: ",start(x)[1],"[",start(x)[2],"] end: ",end(x)[1],"[",end(x)[2],"] freq:",frequency(x)," duration: ",length(x), " samples\r\n")
  for(i in 0:(floor(length(x)/frequency(x))-1) ){
    cat(start(x)[1]+i, "\t", paste(formatC(x[(i*frequency(x)+1): (i*frequency(x)+frequency(x))],width=3,flag="0"), sep="\t\t" ),"\r\n")
  }
  if((length(x)/frequency(x))-floor(length(x)/frequency(x)) != 0){
    i=i+1
    cat(start(x)[1]+i, "\t", paste(formatC(x[(i*frequency(x)+1): length(x)],width=2,flag="0"), sep="\t\t" ),"\r\n")
  }
  invisible(x)
}

cloneDyadStream = function(x, stream){
  if(!"ts"%in%class(x) & !"DyadStream" %in% class(stream)) stop ("a ts object must be cloned with a DyadStream object")
  attributes(x) = c(attributes(x)["tsp"],attributes(stream)[!names(attributes(stream))%in% "tsp"])
  x
}

window.DyadStream = function(x, ...){
  cloneDyadStream(stats::window(x,...),x)
}




 