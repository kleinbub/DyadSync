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
# -patient/therapist are weird for non-clinical dyads
# -generic functions (such as filters, CCFs, etc), should automatically understand the kind of object(= exp, dyad, stream), and behave accordingly
# -experiments shoud be a collection of therapies -> sessions -> signals -> streams (so "therapies" must be added)
# o in alternativa usare dyadId e sessionId per fare sorting e inferire la longitudinalità
#
# REMEMBER:
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
# bestCCF e bestLAG salvati come stream.

DyadExperiment = function(name, sessionList){
  exp =  sessionList
  attributes(exp) = c(attributes(exp), list(name   = name))
  class(exp) = append(class(exp),"DyadExperiment")
  return(exp)
}
is.DyadExperiment = function(x) inherits(x,"DyadExperiment") && length(x)


###
## Nota su dyadsession:
# ci sono 2 possibili disegni:
# a) tante sedute longitudinali di una o più diadi
# b) una seduta di più diadi in varie condizioni/gruppo

DyadSession = function(sessionId,dyadId, signalList=NULL,catList=NULL,start=0){
  ses = list(
    signals = signalList,
    categ = catList
  )
  if(!is.null(signalList)) names(ses$signals) = lapply(ses$signals, function(x) x$name)
  attributes(ses) = c(attributes(ses),list(
    sessionId = sessionId,
    dyadId = dyadId,
    name = paste(sessionId,dyadId,sep="_"),
    start = start
  ))
  class(ses) = append(class(ses),"DyadSession")
  return(ses)
}
is.DyadSession = function(x) inherits(x,"DyadSession") && length(x)


getSession = function(dyadSession){
  attr(dyadSession,"sessionId")
}
getDyad = function(dyadSession){
  attr(dyadSession,"dyadId")
} 
getName = function(dyadSession){
  attr(dyadSession,"name")
} 



DyadCCF = function(#bestCCF, bestLag, 
                   ccfmat, sampRate,
                   settings = list("lagSec"=NULL,"incSec"=NULL,"winSec"=NULL,"accelSec"=NULL,"weight"=NULL,"interpolated"=NULL)){
  CCF =  list(
    
              # bestCCF=bestCCF,
              # bestLag = bestLag,
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
  class(CCF) = append(class(CCF),"DyadCCF")
  return(CCF)
  
}

DyadCategory = function(name, data, sampRate){
  categ = data
  attributes(categ) = c(attributes(categ),list(name = name))
  class(categ) = append(class(categ),"DyadCategory")
  return(categ)
}
is.DyadCategory = function(x) inherits(x,"DyadCategory") && length(x)


DyadSignal = function(name="some signal",patient=NULL,clinician=NULL,sampRate=NULL,dyadNames = c("patient","clinician")){
  sig = list(
    name = name,
    patient =   DyadStream(stream=patient,   name=dyadNames[1], sampRate = sampRate, col = "deeppink3",  lty=1, lwd=2),
    clinician = DyadStream(stream=clinician, name=dyadNames[2], sampRate = sampRate, col = "dodgerblue3", lty=1, lwd=2),
    sampRate = sampRate,
    dyadNames = dyadNames,
    ccf = NULL
  )
  class(sig) = append(class(sig),"DyadSignal")
  return(sig)
} 

is.DyadSignal = function(x) inherits(x,"DyadSignal") && length(x)




### Dyadstream extends ts() including a name and some graphical parameters.
DyadStream = function(stream, name, sampRate=frequency(stream), start=0, col=1, lty=1, lwd=1){
  if(!is(stream,"ts")){
    stream = ts(stream, start = start, frequency = sampRate)
    warning(paste0("Stream '",name,"' was coerced to ts starting at 0, with sampRate of: ",sampRate,"Hz, and duration of: ",end(stream)[1]," seconds."), call.=F)
  }
  attributes(stream) = c(attributes(stream),
                         list(
                            name = name,
                            col=col,
                            lty=lty,
                            lwd = lwd
                        ))
  class(stream) = append("DyadStream",class(stream))
  return(stream)
}

CCFStream = function(stream, lagSec, incSec, winSec, accelSec, weight, interpolated ){
  if(!is.DyadStream(stream))stop("only dyadstream get butter")
  
  attributes(stream) = c(attributes(stream),
                         list(
                            lagSec   = lagSec,
                            incSec   = incSec,
                            winSec   = winSec,
                            accelSec = accelSec,
                            weight   = weight,
                            interpolated = interpolated
                          ))
  class(stream) = append("CCFStream",class(stream))
  return(stream)
}

is.DyadStream = function(x) inherits(x,"DyadStream") && length(x)
is.CCFStream = function(x) inherits(x,"CCFStream") && length(x)

print.DyadStream = function (x, ...) 
{
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

getDyadStreamAttr = function(x){
  a = list(name=attr(x,"name"),col=attr(x,"col"),lty=attr(x,"lty"),lwd=attr(x,"lwd"))
  if(is.CCFStream(x)){a = c(a,list(
    lagSec =attr(x, "lagSec"), incSec =attr(x, "incSec"), winSec =attr(x, "winSec"), accelSec =attr(x, "accelSec"), weight =attr(x, "weight"), interpolated =attr(x, "interpolated")
  ))}
  a
  
}
cloneDyadStream = function(x, stream){
  if(!"ts"%in%class(x) & !"DyadStream" %in% class(stream)) stop ("a ts object must be cloned with a DyadStream object")
  attributes(x) = c(attributes(x),getDyadStreamAttr(stream))
  class(x) = append("DyadStream",class(x))
  if(is.CCFStream(stream)) class(x) = append("CCFStream",class(x))
  x
}

window.DyadStream = function(x, ...){
  cloneDyadStream(stats::window(x,...),x)
}

#generic functions

name <- function(x) {
  attr(x,"name")
}


 