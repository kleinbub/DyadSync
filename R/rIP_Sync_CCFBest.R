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



#' improved ccfBest function with MA and slope filters, to mimic Marci et al 2007 procedure.
#'
#' @param x
#' @param signals
#' @param lagSec
#' @param winSec
#' @param incSec
#' @param accelSec
#' @param slopes      vector of length 2, containing the window size and increments of a "average slope" filter, in seconds
#' @param MA          vector of length 2, containing the window size and increments of the moving average filter, in seconds
#' @param weight_type
#' @param simplify
#' @param outputName
#' @export
ccfBest = function(x, signals="all", lagSec,winSec,incSec,accelSec, slopes=c(NA,NA), MA=c(NA,NA), weight_type=c("center","free","off"),simplify = T, outputName = "CCFBest")
{
  UseMethod("ccfBest",x)
}

#' @export
ccfBest.DyadExperiment = function(experiment, signals="all", lagSec,winSec,incSec,accelSec, slopes=c(NA,NA), MA=c(NA,NA), weight_type=c("center","free","off"),simplify = T, outputName = "CCFBest"){
  cat("\r\nccfBest routine | moving windows cross-correlation with sync maximizing\r\nWith the chosen settings,
      the best lag is able to change by +/-",
      (accelSec)/(incSec),
      "seconds each second of signal." )
  cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]), " follows the ",
             s1Name(experiment[[1]]),"\r\n"))
  if(length(MA)!=2 || (!all(is.numeric(MA)) || !all(is.na(MA)) ) ) MA = c(NA,NA)
  if(winSec%%2==1) warning("Using odd 'winSec' values will cause slightly uncentered correlation windows, due to ts limitations")
  nSessions = length(experiment)
  experiment2 = Map(function(session,iSession){
    if(signals=="all") signals = names(session)
    cat("\r\n",paste(dyadId(session),session(session)))
    session[signals] = Map(function(signal,iSignal){
      cat(" |",name(signal))
      signal = ccfBest(signal, lagSec,winSec,incSec,accelSec, slopes, MA, weight_type, simplify, outputName)
      return(signal)
    }, session[signals])
    #prog(iSession,nSessions)
    return(session)
  },experiment,seq_along(experiment))
  cat("\r\nDone ;)")
  attributes(experiment2)=attributes(experiment)
  return(experiment2)
  if(lagSec > 5)  message("SC latency from stimuli is between 1 and 5 sec. Bigger lags are robably wrong in a stimulus-response perspective")
  
}

#' @export
ccfBest.DyadSignal = function(signal, lagSec,winSec,incSec,accelSec, slopes=c(NA,NA), MA=c(NA,NA), weight_type=c("center","free","off"),simplify = T, outputName = "CCFBest"){
  signal[[outputName]] = CCFBest(NULL,NULL,NULL, lagSec,winSec,incSec,accelSec,weight_type, MA[1], MA[2],slopes)
  origs1 = signal$s1 # save original signals befor applying MA and slope filters
  origs2 = signal$s2
  if(length(slopes)==2 && all(is.numeric(slopes)))
    signal = signalFilter(signal,movAvSlope,slopes[1],slopes[2])
  if(length(MA)==2 && all(is.numeric(MA)))
    signal = signalFilter(signal,movAv,MA[1],MA[2])
  
  signal = dyadCCF(signal,lagSec,winSec,incSec,simplify,outputName)
  signal = vectorBestLag(signal,accelSec,weight_type,outputName)
  signal$s1 = origs1
  signal$s2 = origs2
  
  signal
}


dyadCCF = function(signal,lagSec,winSec,incSec, simplify, outputName = "CCFBest"){
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
  signal[[outputName]]$table = data.frame(matrix(unlist(lcc),ncol=length(ran), byrow = T, dimnames=list(paste0("w",seq_len(n_win)),paste0("lag",ran))))
  colnames(signal[[outputName]]$table) = paste0("lag",ran)
  # signal$ccf$ccfmat[is.na(signal$ccf$ccfmat)] = 0
  xStart = c(start(signal)[1] +trunc(winSec/2),1)
  cat("\r\n freq = ",1/incSec)
  signal[[outputName]]$zero =DyadStream(signal[[outputName]]$table[["lag0"]],"lagZeroSync","#A11F12",start=xStart,frequency=1/incSec )
  #dimostrazione che la CCF inizia a metà della finestra!
  # plot(rangeRescale(window(signal$s1,start=277,end=277+15),-1,1),type="l");abline(v=277+7.5,lty=3)
  # points(signal[[outputName]]$zero)
  # abline(v=277+7.5+incSec*0:5,lty=3)
  return(signal)
}


#vectorBestLag calcola:
#bestLag con i valori in secondi della lag che massimizza la correlazione
#bestCCF con i valori di correlazione a bestLag
#la flessibilità della funzione dipende da:
## -accelSec: indica di quanti secondi può variare il lag ogni secondo di segnale
## - weight_type: BOOL indica se rendere elastica la lag o no (attualmente l'elastico tira verso il valore attuale, non lag0)
vectorBestLag = function(signal, accelSec, weight_type=c("off","center","free"), outputName = "CCFBest"){
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  weight_type = match.arg(weight_type)
  mat = as.matrix(signal[[outputName]]$table)
  sampRate = sampRate(signal)
  
  #ora a è più corretto e funziona anche se simplify = TRUE
  colFrequency = (ncol(mat)-1)/(2*attr(signal[[outputName]],"lagSec")) #quante colonne compongono un secondo?
  incSec = attr(signal[[outputName]],"incSec") #quante righe compongono un secondo?
  #muovere la lag di 1 secondo ogni secondo di segnale, implica muovere di 'acc' colonne per ogni riga:
  acc = colFrequency*incSec
  #a = n colonne per muovere la lag di accelSec secondi al secondo
  a = accelSec*acc
  if(a<1) stop ("With chosen [incSec, sampRate, and simplify], accelSec must be at least:", 1/acc)
  if(a%%1!=0) { #if not integer
    a = max(1,round(a))
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
  signal[[outputName]]$table = cbind(signal[[outputName]]$table, "bestCCF" = bCC, "bestLag" = blag)
  xStart = c(start(signal)[1] +attr(signal[[outputName]],"winSec")/2,1)
  signal[[outputName]]$sync = DyadStream(bCC,  name="bestCCF", "#FF0041", start=xStart, frequency= 1/incSec )
  signal[[outputName]]$lag  = DyadStream(blag, name="bestLAG", "#BFBF00", start=xStart, frequency= 1/incSec )

  return(signal)
}


#marci index ottimizzato. (taglia i valori infiniti a +- 10)
#per far na roba fatta ben si dovrebbe vedere la distribuzione casuale di marciIndex e vedere come e dove tagliare

#' Marci's concordance index
#'
#' @param a a vector of cross-correlations, or a list of them
#' @details see Marci et al. 2007
#' @export
#'

marciIndex = function(a){
  if(is.list(a)){ ##se 'a' è una lista, applica l'indice a ogni elemento della lista
    sapply(a, function(iFile){apply(iFile,2,marciIndex)})
  } else {
    a = na.omit(as.numeric(a))
    x=log(sum(a[a >0]) / abs(sum(a[a<0]))) #formula dell'indice
    if(!is.na(x)){
      if(x>=10)x=10 else if(x <= -10)x=-10
    }
    x
  }
}

#' kleinbub session index
#' A simple summarizing index of synchronization in a given session
#' @details the index calculates the median between high values and multiplies it
#' by the ratio of high values in the time-series.
#' By default, high values are values greater than 0.
#'
#' @param a 
#'
#' @return
#' @export

kleinbubIndex=function(a, threshold = 0){ #KSI = kleinbub session index
  if(is.list(a)){ ##se 'a' è una lista, applica l'indice a ogni elemento della lista
    sapply(a, function(iFile){apply(iFile,2,kleinbubIndex)})
  } else {
  a = na.omit(as.numeric(a))
  b = a[a > threshold]
  ratio = length(b)/length(a)
  median(b) * ratio
  }
}


