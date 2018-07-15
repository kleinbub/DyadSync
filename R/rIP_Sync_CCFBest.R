##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  
# moving windows cross correlation functions
#
############################################################################################################ 
## changelog

############################################################################################################ 
## ToDo

############################################################################################################ 

#calcola la CCF a finestre mobili per diverse lag su cui calcola la CCF elastica
#lagSec = max amount of lag in seconds to analyze for each direction (es: '3' -> -3, -2, -1, 0, 1, 2, 3 )
#winSec, incSec = size and increments of running windows, in seconds.
#accelSec= allowed lag variation (in seconds) between two seconds of the signal
#weight = should the elastic be weighted? Accepts values of 'center','free',or FALSE
#interpolate = if TRUE, the CCF is translated to the same sampling rate of the signal
#simplify = if TRUE, only the ccf matrix will be reduced to discrete seconds of lag.

#' Title
#'
#' @param experiment 
#' @param signals 
#' @param lagSec 
#' @param winSec 
#' @param incSec 
#' @param accelSec 
#' @param weight 
#' @param interpolate 
#' @param simplify 
#'
#' @return
#' @export
#'
#' @examples
#' 

ccfBest = function(x, signals="all", lagSec,winSec,incSec,accelSec,weight_type=c("center","free"),simplify = T)
{
  UseMethod("ccfBest",x)
}

#' @export
ccfBest.DyadExperiment = function(experiment, signals, lagSec,winSec,incSec,accelSec,weight_type,simplify){
  cat("\r\nccfBest routine | moving windows cross-correlation with sync maximizing\r\nWith the chosen settings,
      the best lag is able to change by +/-",
      (accelSec)/(incSec),
      "seconds each second of signal." )
  cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]), " follows the ",
             s1Name(experiment[[1]]),"\r\n"))
  nSessions = length(experiment)
  experiment2 = Map(function(session,iSession){
    if(signals=="all") signals = names(session)
    cat("\r\n",paste(dyadId(session),session(session)))
    session[signals] = Map(function(signal,iSignal){
      cat(" |",name(signal))
      signal = ccfBest(signal, lagSec,winSec,incSec,accelSec,weight_type,simplify)
      return(signal)
    }, session[signals])
    #prog(iSession,nSessions)
    return(session)
  },experiment,seq_along(experiment))
  cat("\r\nDone ;)")
  attributes(experiment2)=attributes(experiment)
  return(experiment2)
}


ccfBest.DyadSignal = function(signal, lagSec,winSec,incSec,accelSec,weight_type,simplify=T){
  if(lagSec > 5)  warning("SC latency from stimuli is between 1 and 5 sec. Bigger lags are robably wrong in a stimulus-response perspective")
  
  signal$CCFBest = CCFBest(NULL,NULL,NULL, lagSec,winSec,incSec,accelSec,weight_type)
  
  signal = dyadCCF(signal,lagSec,winSec,incSec,simplify)
  signal = vectorBestLag(signal,accelSec,weight_type)
  # if(simplify)
  #   signal = reduceCcfLags(signal)
  # if(interpolate){
  #   signal$ccf$ccfmat = winInter(list(signal$ccf$ccfmat), winSec = signal$ccf$settings$winSec, incSec = signal$ccf$settings$incSec, sampRate = signal$sampRate)[[1]]
  #   colnames(signal$ccf$ccfmat) = gsub("[.]","-",colnames(signal$ccf$ccfmat))
  #   signal$ccf$sampRate = signal$sampRate
  #   }
  signal
}


dyadCCF = function(signal,lagSec,winSec,incSec, simplify){
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  sampRate = sampRate(signal)
  #import C correlation function
  C_cor=get("C_cor", asNamespace("stats"))
  
  win = winSec*sampRate
  lagSamp = lagSec*sampRate
  if(!simplify)
    ran = ((-lagSamp)): ((lagSamp))
  else
    ran = seq(-lagSamp,lagSamp, sampRate)
  inc = incSec * sampRate
  iFile = data.frame(signal$s1,signal$s2)
  n_win = ceiling((nrow(iFile)-win-lagSamp+1)/inc) #calcola il numero di finestre
  lcc=lapply(seq_len(n_win)-1, function(iWin)
  { #-----per ciascuna finestra--------
    ab = (iWin*inc +1):(iWin*inc +win+lagSamp) #calcola il range di sample di ciascuna finestra
    xWin = iFile[ab,1];yWin = iFile[ab,2] #estrai i dati della finestra in un vettore per sogg x e y
    # plot(rangeRescale(xWin,1,0),type="l");lines(rangeRescale(yWin,1,0))
    if(max(sum(is.na(xWin)),sum(is.na(yWin))) > length(xWin)/2){ #se ci sono troppi NA, restituisci NA per tutti i lag
      rep(NA, length(ran))
    } else {
      unlist(lapply(ran, function(iLag)
      { #-------per ciascun lag---------
        LAG = abs(iLag)
        xRange = 1:(win);yRange = (LAG+1) :(win+LAG) #applica il lag spostando indietro il secondo soggetto.
        if(iLag<0){k=xRange;xRange=yRange;yRange=k;} #valori alti a lag positivi implicano che il sogg 2 segue sogg 1
        x = xWin[xRange];y = yWin[yRange]
        #         # questi due check permettono di avere un valore quando una finestra ha sempre lo stesso numero
        #         if(var(x,na.rm=T)==0 ) {if(!is.na(x[1])) x[1] = x[1]+ 0.000001 else x[1] = 0.000001}
        #         if(var(y,na.rm=T)==0 ) {if(!is.na(y[1])) y[1] = y[1]+ 0.000001 else y[1] = 0.000001}
        suppressWarnings(cor(x,y,method="pearson",use="complete.obs"))
        if(sum(x!=y,na.rm = T)<2 )  NA #controlla che ci siano almeno 3 punti !=0
        else suppressWarnings(.Call(C_cor, x, y, 2L, FALSE))
      }))
    }#fine else
  }) #fine lapply finestre
  signal$CCFBest$table = data.frame(matrix(unlist(lcc),ncol=length(ran), byrow = T, dimnames=list(paste0("w",seq_len(n_win)),paste0("lag",ran))))
  colnames(signal$CCFBest$table) = paste0("lag",ran)
  # signal$ccf$ccfmat[is.na(signal$ccf$ccfmat)] = 0
  xStart = c(start(signal)[1] +winSec/2,1)
  signal$CCFBest$zero =DyadStream(signal$CCFBest$table[["lag0"]],"lagZeroSync","#A11F12",start=xStart,frequency=1/incSec )
  #dimostrazione che la CCF inizia a metà della finestra!
  # plot(rangeRescale(window(signal$s1,start=277,end=277+15),-1,1),type="l");abline(v=277+7.5,lty=3)
  # points(signal$CCFBest$zero)
  # abline(v=277+7.5+incSec*0:5,lty=3)
  return(signal)
}


#vectorBestLag calcola:
#bestLag con i valori in secondi della lag che massimizza la correlazione
#bestCCF con i valori di correlazione a bestLag
#la flessibilità della funzione dipende da:
## -accelSec: indica di quanti secondi può variare il lag ogni secondo di segnale
## - weight_type: BOOL indica se rendere elastica la lag o no (attualmente l'elastico tira verso il valore attuale, non lag0)
vectorBestLag = function(signal, accelSec, weight_type=c("off","center","free")){
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  weight_type = match.arg(weight_type)
  mat = as.matrix(signal$CCFBest$table)
  sampRate = sampRate(signal)
  
  #ora a è più corretto e funziona anche se simplify = TRUE
  colFrequency = (ncol(mat)-1)/(2*attr(signal$CCFBest,"lagSec")) #quante colonne compongono un secondo?
  incSec = attr(signal$CCFBest,"incSec") #quante righe compongono un secondo?
  #muovere la lag di 1 secondo ogni secondo di segnale, implica muovere di 'acc' colonne per ogni riga:
  acc = colFrequency*incSec
  #a = n colonne per muovere la lag di accelSec secondi al secondo
  a = accelSec*acc
  if(a<1) stop ("With chosen [incSec, sampRate, and simplify], accelSec must be at least:", 1/acc)
  if(a%%1!=0) {
    a = max(1,round(a))
    realAccelSec = a/acc
    #show the warning only once, per setting.
    if(getOption("accelSec_check")!=paste0(accelSec,incSec,sampRate,simplify) ) {
      warning ("With chosen [incSec, sampRate, and simplify], accelSec should be a multiple of ",acc,
               " that returns an integer. Eg:\r\n",paste(1:10/acc,collapse = ", "),
               "\r\nAccelSec was approximated to ",realAccelSec)
    }
    options("accelSec_check"=paste0(accelSec,incSec,sampRate,simplify)) 
  }

  
  lag0 = ceiling(ncol(mat)/2) #la colonna di lag 0
  
  ###DETERMINE WEIGHTS
  if(weight_type!="off"){
    if(weight_type=="free"){
      #wx copre il range possibile dato accSec
      wx = -a:a
    } else if(weight_type=="center"){
      #wx copre tutti i possibili lag, così il peso è sempre centrato su lag0
      cols=ncol(mat)-lag0
      wx = -cols:cols
    }
    weightMalus = 10
    malus = weightMalus  / 100
    maxMalus = 1-weightMalus/100
    weights_val = rangeRescale(dnorm(wx, mean=0, sd=1000),0,malus ) + maxMalus
  }   

  
  first = ifelse(is.na(mat[1,lag0]),NA,lag0) #comincia con lag0 o NA, se mancano le correlazioni nella prima riga
  blag = c(first, rep(0,nrow(mat)-1)) 
  bCC  = c(as.numeric(mat[1,lag0]), rep(0,nrow(mat)-1)) #vettore vuoto comincia col valore di lag0
  
  for(i in 2:nrow(mat)){ #dalla seconda riga fino alla fine della matrice di lag
    if(is.na(mat[i,lag0])){ #if ccf is NA give NA
      bCC[i] = NA
      blag[i] = NA
    } else {
      previous = ifelse(is.na(blag[i-1]), lag0, blag[i-1])
      ran = seq(previous-a,previous+a) #whole range of possible lag change
      realran = ran[ran>0 & ran<= ncol(mat)] #keep only existing columns (>0 and < than ncol)
      possibleVal = mat[i,realran] #the ccf values at 'realran' lags
      
      if (weight_type!="off"){
        if(weight_type=="free"){
          ranweights = weights_val[which(ran>0 & ran<= ncol(mat))]
        } else if(weight_type=="center") {
          ranweights = weights_val[realran] #the weight values of 'realran' lags
        }
        # if(length(which.max(mat[i,realran]*ranweights)+realran[1]-1)==0){ #this is to debug errors
        #   cat("\r\n", mat[i,realran]*ranweights)
        # }
        # if(all(is.na(possibleVal)))possibleVal = rep(1,length(possibleVal)) #this is to debug errors
        
        blag[i] = which.max(mat[i,realran]*ranweights)+realran[1]-1 #the column with highest (weighted) ccf among realran lags
      } else {
        # if(all(is.na(possibleVal))) blag[i] = round(length(possibleVal)/2)
        # else 
        blag[i] = which.max(mat[i,realran])+realran[1]-1
      }
      
      bCC[i]  = mat[i,blag[i]]
    }
    # ##debug tools
    # cat0("\r\nw",lead0(floor((i+15)/60)),":",lead0( (i+15)- floor((i+15)/60)*60) )
    # cat0(" | among ",lead0(realran[1]),"-",lead0(realran[length(realran)]))
    # cat0(" selected: ",lead0(blag[i]), "(",round(bCC[i],3)," t:",(blag[i]-lag0)/signal$sampRate,")")
  }
  
  blag = (blag-lag0)/colFrequency #trasforms blag from columns to seconds
  signal$CCFBest$table = cbind(signal$CCFBest$table, "bestCCF" = bCC, "bestLag" = blag)
  xStart = c(start(signal)[1] +attr(signal$CCFBest,"winSec")/2,1)
  signal$CCFBest$sync = DyadStream(bCC,  name="bestCCF", "#FF0041", start=xStart, frequency= 1/incSec )
  signal$CCFBest$lag  = DyadStream(blag, name="bestLAG", "#BFBF00", start=xStart, frequency= 1/incSec )

  return(signal)
}


#marci index ottimizzato. (taglia i valori infiniti a +- 10)
#per far na roba fatta ben si dovrebbe vedere la distribuzione casuale di marciIndex e vedere come e dove tagliare
#' Marci's concordance index
#'
#' @param a a vector of cross-correlations, or a list of them
#' @details see Marci et al. (2007)
#' @export
#'

marciIndex = function(a){
  if(is.list(a)){ ##se 'a' è una lista, applica l'indice a ogni elemento della lista
    matrix(unlist(lapply(a, function(iFile){apply(iFile,2,marciIndex)})),nrow=length(a),byrow=T)
  } else {
    x=log(sum(a[a >0]) / abs(sum(a[a<0]))) #formula dell'indice
    if(!is.na(x)){
      if(x>=10)x=10 else if(x <= -10)x=-10
    }
    x
  }
}

#this function extracts a column from ccfmat of a given signal
#the data is extracted as a DyadStream object with appropriate frequency
#' Title
#'
#' @param signal 
#' @param lag 
#' @param col 
#' @param lty 
#' @param lwd 
#' @param interpolate 
#'
#' @return
#' @export
#'
#' @examples
#' 
# getCCF = function(signal, lag, col="red", lty=2,lwd=2, interpolate = F){
#   stop("  *getCCF è obsoleta tranne per estrarre dati dalla matrice (es: lag0). e va aggiornata.")
#   if(!is.DyadSignal(signal)) stop("Only objects of class DyadSignal can be processed by this function")
#   # if(!grepl("lag",lag, ignore.case = T) & grepl("best",lag, ignore.case = T)) {
#   #   lag = "bestCCF"
#   # }  else if (grepl("lag",lag, ignore.case = T) & grepl("best",lag, ignore.case = T)) {
#   #   lag = "bestLag"
#   # } else if(suppressWarnings(!is.na(as.numeric(lag)))){
#   #   lag = paste0("lag",as.character(lag))
#   # }
#   if(suppressWarnings(!is.na(as.numeric(lag))))
#     lag = paste0("lag",as.character(lag))
#   if (! lag %in% colnames(signal$ccf$ccfmat)) stop("Specified lag value '",lag,"' was not found. Available lags: ", paste(colnames(signal$ccf$ccfmat), collapse=" | "))
#   #cat0("\r\n CCF column '",lag,"' was selected")
#   if(!signal$ccf$settings$interpolated){
#     if(interpolate){
#       res = winInter(signal$ccf$ccfmat[as.character(lag)],
#                      winSec = signal$ccf$settings$winSec,
#                      incSec = signal$ccf$settings$incSec,
#                      sampRate = signal$sampRate)
#       res = unlist(res)
#       newInterpolate = T
#       signal$ccf$sampRate = signal$sampRate
#     } else {
#       warning("CCF is not interpolated. Sample rate setting might be wrong, or manifest other unexpected result",call.=F)
#       res = signal$ccf$ccfmat[as.character(lag)]
#       newInterpolate = F
#     } 
#     
#   } else {
#     res = signal$ccf$ccfmat[as.character(lag)]
#     newInterpolate = T
#   }
#   settings = list ( 
#     lagSec       = signal$ccf$settings$lagSec,
#     incSec       = signal$ccf$settings$incSec,
#     winSec       = signal$ccf$settings$winSec,
#     accelSec     = signal$ccf$settings$accelSec,
#     weight       = signal$ccf$settings$weight,
#     interpolated = newInterpolate)
#   stream = DyadStream(ts(res, frequency = signal$ccf$sampRate,start=0),name = paste0("CCF - ",lag), settings= settings, col = col, lty = lty, lwd = lwd)
#   
# }


### STOP! don't calculate them if you don't need them! ####
# reduceCcfLags = function(signal){
#   if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
#   if(signal$sampRate > 1){
#     best = signal$ccf$ccfmat[,c("bestCCF","bestLag")]
#     ran = seq(-signal$sampRate*signal$ccf$settings$lagSec,signal$sampRate*signal$ccf$settings$lagSec,signal$sampRate)
#     ran = ran + signal$sampRate*signal$ccf$settings$lagSec + 1
#     signal$ccf$ccfmat = cbind(signal$ccf$ccfmat[,ran], best)
#     return(signal)
#   } else {
#     return(signal)
#   }
# }