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
#accelSec= allowed lag variation (in seconds) between two seconds of the x
#weight = should the elastic be weighted? Accepts values of 'center','free',or FALSE
#interpolate = if TRUE, the CCF is translated to the same sampling rate of the x
#simplify = if TRUE, only the ccf matrix will be reduced to discrete seconds of lag.



#' improved ccfBest function with MA and slope filters, to mimic Marci et al 2007 procedure.
#'
#' @param x
#' @param signal
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
ccfBest = function(x, signal="SC", lagSec,winSec,incSec,accelSec, slopes=c(NA,NA), MA=c(NA,NA), weight_type=c("center","free","off"),simplify = T, outputName = "CCFBest")
{
  UseMethod("ccfBest",x)
}

#' @export
ccfBest.DyadExperiment = function(experiment, signal="all", lagSec,winSec,incSec,accelSec, slopes=c(NA,NA), MA=c(NA,NA), weight_type=c("center","free","off"),simplify = T, outputName = "CCFBest"){

  if(length(MA)!=2 || (!all(is.numeric(MA)) || !all(is.na(MA)) ) ) MA = c(NA,NA)
  # if(winSec%%2==1) warning("Using odd 'winSec' values will cause slightly uncentered correlation windows, due to ts limitations")
  nSessions = length(experiment)
  
  ###PARALLELIZE
  
  # progresbar
  nSessions = length(experiment)
  pb <- progress::progress_bar$new(
    format = "Calculation::percent [:bar] :elapsed | ETA: :eta",
    total = nSessions,    # number of iterations
    width = 60, 
    show_after=0 #show immediately
  )
  progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  } 
  opts <- list(progress = progress)
  
  warningsFile = "CCFWarnings"
  if (file.exists(warningsFile)) {
    unlink(warningsFile)
  }
  cores=parallel::detectCores()-1
  cat(paste0("\r\nPerforming parallelized computation of ccfBbest using ",cores," cores.\r\n")) #verified!
  cat("With the chosen settings,
      the best lag is able to change by +/-",
      (accelSec)/(incSec),
      "seconds each second of x." )
  cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]), " follows the ",
             s1Name(experiment[[1]]),"\r\n"))
  cl <- parallel::makeCluster(cores[1], outfile=warningsFile) #not to overload your computer
  doSNOW::registerDoSNOW(cl)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  pb$tick(0)
  experiment2 <- foreach::foreach(
    iSession = 1:nSessions,
    .options.snow = opts, .errorhandling='pass'
  ) %dopar% { 
    xsignal = experiment[[iSession]][[signal]]
    xsignal = ccfBest(xsignal, lagSec,winSec,incSec,accelSec, slopes, MA, weight_type, simplify, outputName)
    experiment[[iSession]][[signal]] = xsignal
    experiment[[iSession]]
    
  }
  parallel::stopCluster(cl) 
  cat("\r\nDone ;)\r\n")
  
  #restore session names and attributes
  for(iSession in 1:nSessions){
    if(!is.null(experiment2[[iSession]]$message)){
      write(paste0("QQQZ",experiment2[[iSession]]$message,"\n"),file=warningsFile,append=TRUE)
    }
    experiment2[[iSession]] = cloneAttr(experiment[[iSession]],experiment2[[iSession]])
  }

  #restore class and attributes
  experiment2 = cloneAttr(from=experiment,to=experiment2)
  names(experiment2) = names(experiment)
  
  if (file.exists(warningsFile)) {
    wr = readLines(warningsFile)
    wr = wr[grepl("QQQZ",wr,fixed = T)]
    if(length(wr)>0){
      cat("\r\nIssues:\r\n")
      for(i in 1:length(wr)){
        cat(substr(wr[i], start = 5, stop=10000), "\r\n")
      }
      stop("Please fix the issues before proceding.")
      
    }
    unlink(warningsFile)
  }
  
  return(experiment2)
  
}

#' @export
ccfBest.DyadSignal = function(x, lagSec,winSec,incSec,accelSec, slopes=c(NA,NA), MA=c(NA,NA), weight_type=c("center","free","off"),simplify = T, outputName = "CCFBest"){
  x[[outputName]] = CCFBest(NULL,NULL,NULL, lagSec,winSec,incSec,accelSec,weight_type, MA[1], MA[2],slopes)
  origs1 = x$s1 # save original signals befor applying MA and slope filters
  origs2 = x$s2
  if(length(slopes)==2 && all(is.numeric(slopes)))
    x = signalFilter(x,movAvSlope,slopes[1],slopes[2])
  if(length(MA)==2 && all(is.numeric(MA)))
    x = signalFilter(x,movAv,MA[1],MA[2])
  
  x = dyadCCF(x,lagSec,winSec,incSec,simplify,outputName)
  x = vectorBestLag(x,accelSec,weight_type,outputName)
  x$s1 = origs1
  x$s2 = origs2
  
  x
}


dyadCCF = function(x,lagSec,winSec,incSec, simplify, outputName = "CCFBest"){
  if(!is(x,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  SR = frequency(x)
  #import C correlation function
  C_cor=get("C_cor", asNamespace("stats"))
  
  win = winSec*SR
  lagSamp = lagSec*SR
  if(!simplify)
    ran = ((-lagSamp)): ((lagSamp))
  else
    ran = seq(-lagSamp,lagSamp, SR)
  inc = incSec * SR
  iFile = data.frame(x$s1,x$s2)
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
  x[[outputName]]$table = data.frame(matrix(unlist(lcc),ncol=length(ran), byrow = T, dimnames=list(paste0("w",seq_len(n_win)),paste0("lag",ran))))
  colnames(x[[outputName]]$table) = paste0("lag",ran)
  # x$ccf$ccfmat[is.na(x$ccf$ccfmat)] = 0
  xStart = c(start(x) +trunc(winSec/2),1)
  x[[outputName]]$zero =rats(x[[outputName]]$table[["lag0"]],start=xStart,frequency=1/incSec, timeUnit="second", unit="Pearson correlation" )
  #dimostrazione che la CCF inizia a metà della finestra!
  # plot(rangeRescale(window(x$s1,start=277,end=277+15),-1,1),type="l");abline(v=277+7.5,lty=3)
  # points(x[[outputName]]$zero)
  # abline(v=277+7.5+incSec*0:5,lty=3)
  return(x)
}


#vectorBestLag calcola:
#bestLag con i valori in secondi della lag che massimizza la correlazione
#bestCCF con i valori di correlazione a bestLag
#la flessibilità della funzione dipende da:
## -accelSec: indica di quanti secondi può variare il lag ogni secondo di segnale
## - weight_type: BOOL indica se rendere elastica la lag o no (attualmente l'elastico tira verso il valore attuale, non lag0)
vectorBestLag = function(x, accelSec, weight_type=c("off","center","free"), outputName = "CCFBest"){
  if(!is(x,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  weight_type = match.arg(weight_type, choices = c("off","center","free"))
  mat = as.matrix(x[[outputName]]$table)
  SR = frequency(x)
  
  #ora a è più corretto e funziona anche se simplify = TRUE
  colFrequency = (ncol(mat)-1)/(2*attr(x[[outputName]],"lagSec")) #quante colonne compongono un secondo?
  incSec = attr(x[[outputName]],"incSec") #quante righe compongono un secondo?
  #muovere la lag di 1 secondo ogni secondo di segnale, implica muovere di 'acc' colonne per ogni riga:
  acc = colFrequency*incSec
  #a = n colonne per muovere la lag di accelSec secondi al secondo
  a = accelSec*acc
  if(a<1) stop ("With chosen [incSec, SR, and simplify], accelSec must be at least:", 1/acc)
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
  }
  
  blag = (blag-lag0)/colFrequency #trasforms blag from columns to seconds
  x[[outputName]]$table = cbind(x[[outputName]]$table, "bestCCF" = bCC, "bestLag" = blag)
  xStart = c(start(x) +attr(x[[outputName]],"winSec")/2,1)
  x[[outputName]]$sync = rats(bCC,  start=xStart, frequency= 1/incSec, timeUnit = "second", unit = "Pearson correlation" )
  x[[outputName]]$lag  = rats(blag, start=xStart, frequency= 1/incSec, timeUnit = "second", unit = "seconds" )

  return(x)
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


