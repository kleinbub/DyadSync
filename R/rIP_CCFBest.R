##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  
## expCCF.R
# moving windows cross correlation functions
#
############################################################################################################ 
## changelog
# v1.31 - peak-centered-flex-CCF v1.93 implemented
# v1.3  - peak-centered-flex-CCF v1.2 implemented
# v1.21 - "sessions" removed. all ts() now have explicit start = 0
# v1.2  - vectorBestLag() pimped:
#           -weight is no more T/F but can get values of 'center' (returning toward lag0)
#            or 'free' (returning to last best value)
#           -bestLag is now exported in seconds, not more in columns
#       - compatible with DyadStream 1.2
#       - added getCCF to export elements of the ccfMat to DyadStream format
# v1.1  - Aggiunti commenti e print ausiliari. 
#
############################################################################################################ 
## ToDo
# - explore MABestCCF with LOESS regression instead of Moving average
# - IMPORTANT: lag calculations should use full windows. Thus dyadCCF() should use larger than
#   set windows and use only adequate points
# - IMPORTANT: bestCCF should ensure biunivocity (should it?)
# - bestCCF could do 2 passages, one forward and one backward and than average or weight the results,
#   to avoid suboptimal first-come-first-served situations
#
############################################################################################################ 

#calcola la CCF a finestre mobili per diverse lag su cui calcola la CCF elastica
#lagSec = max amount of lag in seconds to analyze for each direction (es: '3' -> -3, -2, -1, 0, 1, 2, 3 )
#winSec, incSec = size and increments of running windows, in seconds.
#accelSec= allowed lag variation (in seconds) between tow consecutive ccf windows (= lag variation / 1 incSec)
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

ccfBest = function(x, signals="all", lagSec,winSec,incSec,accelSec,weight=c("center","free"),interpolate=T,simplify = T)
{
  UseMethod("ccfBest",x)
}

#' @export
ccfBest.DyadExperiment = function(experiment, signals, lagSec,winSec,incSec,accelSec,weight,interpolate,simplify){
  xSampRate  = sampRate(experiment[[1]][[ifelse(signals=="all",1,signals)]])
  cat("\r\nccfBest routine | moving windows cross-correlation with sync maximizing\r\nWith the chosen settings,
      the best lag is able to change by +/-",
      (accelSec*xSampRate)/(incSec*xsampRate),
      "seconds each second of signal." )
  cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]), " follows the ",
             s1Name(experiment[[1]]),"\r\n"))
  nSessions = length(experiment)
  experiment2 = Map(function(session,iSession){
    if(signals=="all") signals = names(session)
    cat("\r\n",paste(id(session),session(session)))
    session[signals] = Map(function(signal,iSignal){
      cat(" |",name(signal))
      signal = ccfBest(signal, lagSec,winSec,incSec,accelSec,weight,interpolate,simplify)
      return(signal)
    }, session[signals])
    #prog(iSession,nSessions)
    return(session)
  },experiment,seq_along(experiment))
  cat("\r\nDone ;)")
  attributes(experiment2)=attributes(experiment)
  return(experiment2)
}


signalCCF = function(signal, lagSec,winSec,incSec,accelSec,weight,interpolate=T,simplify=T){
  ## debug
  # lagSec = 5
  # winSec = 10
  # incSec = 1
  # accelSec = 0.3
  # weight= F
  # interpolate=T
  # simplify=F
  #####
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")

  #costruttore di classe
  signal$ccf = CCFMatrix(NULL, sampRate = incSec, lagSec=lagSec,incSec=incSec,winSec=winSec,accelSec=accelSec)
  #funzione vera e propria
    signal = dyadCCF(signal)
    signal = vectorBestLag(signal)
    # signal = ppBestLag(signal)

  
  # signal = MABestLag(signal);warning("experimental feature, get back to vectorBestLag")
  if(simplify)
    signal = reduceCcfLags(signal)
  if(interpolate){
    signal$ccf$ccfmat = winInter(list(signal$ccf$ccfmat), winSec = signal$ccf$settings$winSec, incSec = signal$ccf$settings$incSec, sampRate = signal$sampRate)[[1]]
    colnames(signal$ccf$ccfmat) = gsub("[.]","-",colnames(signal$ccf$ccfmat))
    signal$ccf$sampRate = signal$sampRate
    }
  signal
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
getCCF = function(signal, lag, col="red", lty=2,lwd=2, interpolate = F){
  if(!is.DyadSignal(signal)) stop("Only objects of class DyadSignal can be processed by this function")
  # if(!grepl("lag",lag, ignore.case = T) & grepl("best",lag, ignore.case = T)) {
  #   lag = "bestCCF"
  # }  else if (grepl("lag",lag, ignore.case = T) & grepl("best",lag, ignore.case = T)) {
  #   lag = "bestLag"
  # } else if(suppressWarnings(!is.na(as.numeric(lag)))){
  #   lag = paste0("lag",as.character(lag))
  # }
  if(suppressWarnings(!is.na(as.numeric(lag))))
      lag = paste0("lag",as.character(lag))
  if (! lag %in% colnames(signal$ccf$ccfmat)) stop("Specified lag value '",lag,"' was not found. Available lags: ", paste(colnames(signal$ccf$ccfmat), collapse=" | "))
  #cat0("\r\n CCF column '",lag,"' was selected")
  if(!signal$ccf$settings$interpolated){
    if(interpolate){
      res = winInter(signal$ccf$ccfmat[as.character(lag)],
                     winSec = signal$ccf$settings$winSec,
                     incSec = signal$ccf$settings$incSec,
                     sampRate = signal$sampRate)
      res = unlist(res)
      newInterpolate = T
      signal$ccf$sampRate = signal$sampRate
    } else {
      warning("CCF is not interpolated. Sample rate setting might be wrong, or manifest other unexpected result",call.=F)
      res = signal$ccf$ccfmat[as.character(lag)]
      newInterpolate = F
    } 
    
  } else {
    res = signal$ccf$ccfmat[as.character(lag)]
    newInterpolate = T
    }
    settings = list ( 
            lagSec       = signal$ccf$settings$lagSec,
            incSec       = signal$ccf$settings$incSec,
            winSec       = signal$ccf$settings$winSec,
            accelSec     = signal$ccf$settings$accelSec,
            weight       = signal$ccf$settings$weight,
            interpolated = newInterpolate)
  stream = DyadStream(ts(res, frequency = signal$ccf$sampRate,start=0),name = paste0("CCF - ",lag), settings= settings, col = col, lty = lty, lwd = lwd)

}



reduceCcfLags = function(signal){
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  if(signal$sampRate > 1){
    best = signal$ccf$ccfmat[,c("bestCCF","bestLag")]
    ran = seq(-signal$sampRate*signal$ccf$settings$lagSec,signal$sampRate*signal$ccf$settings$lagSec,signal$sampRate)
    ran = ran + signal$sampRate*signal$ccf$settings$lagSec + 1
    signal$ccf$ccfmat = cbind(signal$ccf$ccfmat[,ran], best)
    return(signal)
  } else {
    return(signal)
  }
}


dyadCCF = function(signal){
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  lagSec = signal$ccf$settings$lagSec
  winSec = signal$ccf$settings$winSec
  incSec = signal$ccf$settings$incSec
  sampRate = signal$sampRate
  
  #cat("\r\nMy Setting: lagSec = ",lagSec,"; winSec:",winSec,"; incSec:",incSec,"; sampRate:",sampRate)
  #cat("\r\ndyad CCF function -- v.1.2\r\ninput:",length(fileList),"dyads.\r\nFrequency:",sampRate)
  #   xname = signal$dyadNames[1];yname = signal$dyadNames[2];
  #   cat(paste0("\r\nHigh Sync at positive lags implies that the ",yname, " follows the ",xname,"\r\n"))
  win = winSec*sampRate
  lagSamp = lagSec*sampRate
  #if(lagUnit == "sample")
  ran = ((-lagSamp)): ((lagSamp))
  #else
  #  ran = seq(-lagSamp,lagSamp, sampRate)
  inc = incSec * sampRate
  #if(winSec!=incSec) warning("With overlapping windows, the sync series has to be shifted forward by half window in order to be aligned to the original signal!\r\n")
  iFile = data.frame(signal$s1,signal$s2)
  
  n_win = ceiling((nrow(iFile)-win-lagSamp+1)/inc) #calcola il numero di finestre per ciascun file
  lcc=lapply(seq_len(n_win)-1, function(iWin)
  { #-----per ciascuna finestra--------
    ab = (iWin*inc +1):(iWin*inc +win+lagSamp) #calcola il range di sample di ciascuna finestra
    xWin = iFile[ab,1];yWin = iFile[ab,2] #estrai i dati della finestra in un vettore per sogg x e y
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
      }))
    }#fine else
  }) #fine lapply finestre
signal$ccf$ccfmat = data.frame(matrix(unlist(lcc),ncol=length(ran), byrow = T, dimnames=list(paste0("w",seq_len(n_win)),paste0("lag",ran))))
colnames(signal$ccf$ccfmat) = paste0("lag",ran)
signal$ccf$ccfmat[is.na(signal$ccf$ccfmat)] = 0
return(signal)
}


#vectroBestLag aggiunge due colonne a ccfmat:
#bestLag con i valori in secondi della lag che massimizza la correlazione
#bestCCF con i valori di correlazione a bestLag
#la flessibilità della funzione dipende dalle impostazioni signal$setting del segnale:
## -accelSec: indica di quanti secondi può variare il lag, fra due valori di CCF.
## - weight: BOOL indica se rendere elastica la lag o no (attualmente l'elastico tira verso il valore attuale, non lag0)
vectorBestLag = function(signal){
  ##NB AccelSec indica di quanti secondi può variare il lag, fra due valori di CCF.
  ## Nota che il tempo che intercorre fra due valori di CCF è incSec.
  ## ToDo: il problema è che com funziona ora l'elastico tende a preservare lo status quo, mentre dovrebbe tirare verso lag0.
  
  #   ##debug
  # qwe = expCCF(lr, "SC", lagSec = 5, winSec = 30, incSec = 1, accelSec = 0.3, weight = "center", interpolate = F,simplify = F)
  # signal = qwe$sessions$seduta01_2015$SC
  # signal$ccf$ccfmat =  signal$ccf$ccfmat[1:(ncol(signal$ccf$ccfmat)-2)]
  # rm(qwe)
  
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  accelSec = signal$ccf$settings$accelSec
  weight_type = signal$ccf$settings$weight
  if(!is.logical(weight_type)) weight = TRUE else if(!weight_type) weight = FALSE else stop("Weight setting not recognized. Please use either 'free' or 'center'.")
  
  mat = as.matrix(signal$ccf$ccfmat)
  lagsPerSec = signal$sampRate
  if(ncol(mat) > 5*2*lagsPerSec+1)  warning("SC latency from stimuli is between 1 and 5 sec. Bigger lags are robably wrong for most signals")
  
  ## a is the amount of columns in which the system is allowed to search each row
  ## each column is 1/sampRate of a second
  ## each row is incSec seconds
  ## so the best lag can variate by accelSec/incSec seconds of lag for each second of raw signal.
  ## or in other terms, by accelSec seconds of lag, each window increment.
  a = round(lagsPerSec * accelSec)

  lag0 = ceiling(ncol(mat)/2) #la colonna di lag 0
  
  ###DETERMINE WEIGHTS
  if(weight){
    if(grepl("free",weight_type,ignore.case = T)){
      #questo applica i pesi centrati su a. forse sarebbe il caso di applicarli centrati su lag0?
      weiz = dnorm(-a:a,mean = 0, sd = a*0.8) + (1-max(dnorm(-a:a,mean = 0, sd = a*0.8)))
    } else if(grepl("center",weight_type, ignore.case = T)){
      #prova:
      cols=ncol(mat)-lag0
      #con queste impostazioni lag +- 5 secondi è penalizzata del 10%
      #indipendentemente dal numero di lag considerati (sd costante)
      weights_val = dnorm(-cols:cols,mean = 0, sd = 30)*10
      weights_val = weights_val + 1-max(weights_val)
      #cols è il numero di lag * sampRate
      #quindi per 5s di lag a 10 Hz = 50
      #plot( -cols:cols,weights_val,ylim=c(0.8,1))
      #plot( -50:50, dnorm(-50:50,mean = 0, sd = 44)*100,ylim=c(0,1))
    } else stop("Weight setting not recognized. Please use either 'free' or 'center'.")
  }

  if(is.na(mat[1,lag0])){
    blag = c(NA, rep(0,nrow(mat)-1)) #vettore vuoto comincia con NA
  } else   blag = c(lag0, rep(0,nrow(mat)-1)) #vettore vuoto comincia a lag0
  bCC  = c(as.numeric(mat[1,lag0]), rep(0,nrow(mat)-1)) #vettore vuoto comincia col valore di lag0
  
  for(i in 2:nrow(mat)){ #dalla seconda riga fino alla fine della matrice di lag

  ##### debug
  # for(i in 475:550){
    # bCC = bCC[475:550]
    # blag = blag[475:550]
    
    if(is.na(mat[i,lag0])){ #if ccf is NA give NA
      bCC[i] = NA
      blag[i] = NA
    } else {
      previous = ifelse (is.na(blag[i-1]), lag0, blag[i-1])
      ran = seq(previous-a,previous+a) #whole range of possible lag change
      realran = ran[ran>0 & ran<= ncol(mat)] #keep only existing columns (>0 and < than ncol)
      possibleVal = mat[i,realran] #the ccf values at 'realran' lags

      if (weight){
        if(grepl("free",weight_type,ignore.case = T)){
          ranweights = weiz[which(ran>0 & ran<= ncol(mat))]
        } else ranweights = weights_val[realran] #the weight values of 'realran' lags
        
        if(length(which.max(mat[i,realran]*ranweights)+realran[1]-1)==0){ #this is to debug errors
          cat("\r\n", mat[i,realran]*ranweights)
        }
        if(all(is.na(possibleVal)))possibleVal = rep(1,length(possibleVal)) #this is to debug errors
        
        blag[i] = which.max(mat[i,realran]*ranweights)+realran[1]-1 #the column with highest (weighted) ccf among realran lags
      } else {
        if(all(is.na(possibleVal))) blag[i] = round(length(possibleVal)/2)
        else blag[i] = which.max(mat[i,realran])+realran[1]-1
      }
        
      bCC[i]  = mat[i,blag[i]]
    }
    # ##debug tools
    # cat0("\r\nw",lead0(floor((i+15)/60)),":",lead0( (i+15)- floor((i+15)/60)*60) )
    # cat0(" | among ",lead0(realran[1]),"-",lead0(realran[length(realran)]))
    # cat0(" selected: ",lead0(blag[i]), "(",round(bCC[i],3)," t:",(blag[i]-lag0)/signal$sampRate,")")
  }
  
  blag = (blag-lag0)/signal$sampRate #trasforms blag from columns to seconds
  signal$ccf$ccfmat = cbind(signal$ccf$ccfmat, "bestCCF" = bCC, "bestLag" = blag)

  
  return(signal)
}


##MABestLag è una funzione sperimentale che invece di seguire un percorso vector, trova prima i valori di lag che massimizzano
## la correlazione, e dopo crea uno smooth tramite finestra mobile
## i risultati non sono molto incoraggianti, eventualmente provare con un sistema di smoothing LOWESS

MABestLag = function(signal){
 
  #   ##debug
  # qwe = expCCF(lr, "SC", lagSec = 5, winSec = 30, incSec = 1, accelSec = 0.3, weight = "center", interpolate = F,simplify = F)
  # signal = qwe$sessions$seduta01_2015$SC
  # signal$ccf$ccfmat =  signal$ccf$ccfmat[1:(ncol(signal$ccf$ccfmat)-2)]
  # rm(qwe)
  
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  mat = as.matrix(signal$ccf$ccfmat)
  lagsPerSec = signal$sampRate
  if(ncol(mat) > 5*2*lagsPerSec+1)  warning("SC latency from stimuli is between 1 and 5 sec. Bigger lags are robably wrong for most signals")
  lag0 = ceiling(ncol(mat)/2) #la colonna di lag 0
  
  ###DETERMINE WEIGHTS
  weights_val = dnorm(-cols:cols,mean = 0, sd = 30)*30
  weights_val = weights_val + 1-max(weights_val)
  #apply weights
  wmat = t(apply(mat,1,function(x) x*weights_val))
  #find best overall lag
  a = apply(wmat,1, which.max) 

  ##ora media mobile del lag, e poi ritrova le colonne giuste su bCC
  win = 30
  win2 = floor(win/2)
  len = length(a)
  f = unlist(lapply(seq_along(a),function(t){
    i1 = ifelse(t <= win2, 1, (t-win2) )
    i2 = ifelse(t+win2 >= len, len, (t+win2) )
    sum(a[i1:i2])/length(i1:i2)
  }))

  blag = round(f)
  bCC = mat[blag]
  blag = (blag-lag0)/signal$sampRate #trasforms blag from columns to seconds
  signal$ccf$ccfmat = cbind(signal$ccf$ccfmat, "bestCCF" = bCC, "bestLag" = blag)
  return(signal)
}