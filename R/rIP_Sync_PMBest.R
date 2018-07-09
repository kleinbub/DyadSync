##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  

#
############################################################################################################ 
## changelog
# v3.0 compatibile (ma non testato) con rIP_classes v3.0
# v2.0 integrato in rIP package
#
############################################################################################################ 
## ToDo

############################################################################################################ 




#' expPPsync
#'
#' @param experiment 
#' @param signals 
#' @param lagSec 
#' @param sgol_p 
#' @param sgol_n 
#' @param weightMalus weightmalus è la percentuale di malus per il lag più estremo. Es, con weightMalus= 20 e r = 1 al massimo lag, la trasformazione diventa r' = 0.8
#'
#' @return
#' @export
#'
#' @examples
expPPsync = function(experiment, signals="all", lagSec=7, sgol_p = 2, sgol_n = 25, weightMalus = 30){
  if(!is(experiment,"DyadExperiment")) stop("Only objects of class DyadExperiment can be processed by this function")
  cat(paste0("\r\nAre you sure about lag direction?\r\n"))
  #   cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]),
  #   " follows the ", s1Name(experiment[[1]]),"\r\n"))
  nSessions = length(experiment)
  experiment2 = Map(function(session,iSession){
    if(signals=="all") signals = names(session)
    cat("\r\n",paste(dyadId(session),session(session)))
    session[signals] = Map(function(signal){
      cat(" |",name(signal))
      #the first function calculates best lag
      signal = ppBest(signal,lagSec = lagSec, sgol_p = sgol_p, sgol_n = sgol_n, weightMalus = weightMalus)
      #the second calculates the actual sync values
      signal = ppSync(signal)
      return(signal)
    }, session[signals])
    return(session)
  },experiment,seq_along(experiment))
  cat("\r\nDone ;)")
  classAttr(experiment2)=classAttr(experiment)
  return(experiment2)
}


## peak picking best lag
ppBest = function(signal,lagSec=8, sgol_p = 2, sgol_n = 25, weightMalus = 30) {
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")

  d = signal$s2
  d2  = signal$s1
  sampRate = sampRate(signal)
  ransamp = lagSec * sampRate
  
  ### peaks-valleys detection ########
  s1p = peakFinder(d,  sgol_p, sgol_n, mode = "p", 0.1)
  s1v = peakFinder(d,  sgol_p, sgol_n, mode = "v", 0.5)
  s2p = peakFinder(d2, sgol_p, sgol_n, mode = "p", 0.1)
  s2v = peakFinder(d2, sgol_p, sgol_n, mode = "v", 0.5)
  allpik1 = sort(c(s1p$samples,s1v$samples))       #the sample positions containing a peak or valley in d
  allpik2 = sort(c(s2p$samples,s2v$samples))    #the same in d2
  allBool2 = (s2p$bool + s2v$bool) != 0 #the same in d2
  allSec1 = time(d)[allpik1]
  allSec2 = time(d2)[allpik2]
  
  #full matrix with all matches values
  M = matrix(NA ,nrow=length(allpik1),ncol=length(allpik2))
  
  for(i in 1:length(allpik1)){
    ipik = allpik1[i] 
    #trova il range attorno al picco i in cui cercare la lag
    ###NB nella v1.1 questo avviene da valle a valle, per la versione fissa guarda v1.0c
    search_range = (ipik-ransamp):(ipik+ransamp)
    ## seleziona il range di confronto va dal picco/valle precedente a quello successivo di ipik
    #se non ci sono picchi/valli prima, parti dall'inizio del segnale
    a = ifelse(any(allpik1<ipik),  max(allpik1[allpik1<ipik]), 1)
    #se non ci sono picchi/valli dopo, usa la fine del segnale
    b = ifelse(any(allpik1>ipik),  min(allpik1[allpik1>ipik]), length(d))
    ab = a:b
    
    search_range[search_range<=0] = NA
    ab[ab<=0] = NA
    matches = which(allBool2[search_range]) # trova picchi in d2 nell'intorno di ipik (su d)
    matches = matches + (ipik-ransamp) -1     # passa da posizione relativa a search_range a posizione assoluta su d2
    matches_i = which(allpik2 %in% matches)   
    nMatch = length(matches)
    
    if(nMatch>0) {
      ###per ogni match
      for(f in 1:length(matches)){
        ipik2 = matches[f]
        lagv = matches[f] - ipik #distanza fra i due picchi, in samples. Valori negativi indicano giallo anticipa blu
        
        #trova il range valle-valle o picco-picco del segno matchato
        a2 = ifelse(any(allpik2<ipik2),  max(allpik2[allpik2<ipik2]), 1)
        b2 = ifelse(any(allpik2>ipik2),  min(allpik2[allpik2>ipik2]), length(d2))
        ab2= a2:b2

        ## MODE:2 
        ## applica la lag, poi tieni solo l'intersezione dei due segni 
        ## d1:     v---p--------v
        ## d2:  v------p----v
        ## keep:   |--------|
        ## stessa lunghezza (laggato) dell'altro picco.
        l_mar = min(ipik2-a2, ipik-a) #distanza dal picco alla valle sx
        r_mar = min(b2-ipik2, b-ipik) #distanza dal picco alla valle dx
        range1 = (ipik -l_mar):(ipik +r_mar)
        range2 = (ipik2-l_mar):(ipik2+r_mar)
        
        # plot(ab,d[ab],type="l",ylim=c(6,7.5),xlim=c(320,400)); lines(ab2 - lagv,d2[ab2])
        # abline(v=ipik); abline(v=ipik2 - lagv);abline(v=ipik-l_mar); abline(v=ipik+r_mar)
        # lines(range1, d[range1], col="red",lty=2,lwd=2);lines(range2 - lagv, d2[range2], col="red",lty=2,lwd=2)
        
        toCor1 = d[range1]
        toCor2 = d2[range2]
        
        #invece delle correlazioni usa la differenza fra le derivate normalizzate
        m2d = 1-mean(abs(rangeRescale(diff(toCor1),0,1) - rangeRescale(diff(toCor2),0,1)))
        thisCor = m2d

        ## WEIGHT CORRELATION
        # weightMalus = 50 #percentuale di riduzione per il lag più estremo
        malus = weightMalus  / 100
        maxMalus = 1-weightMalus/100
        wx = (-ransamp):+ransamp
        weights_val = rangeRescale(dnorm(wx, mean=0, sd=1000),0,malus ) + maxMalus
        weightCor = thisCor * weights_val[lagv+ransamp+1]
        
        ## POPULATE FINAL MATRIX
        M[i,matches_i[f]] = weightCor
      }
    }
    
  }
  
  ## SORTING APPROACH n2: try all combinations
  ## https://cs.stackexchange.com/questions/91502
  M[is.na(M)] = 0
  ## ignora le correlazioni negative. Non sono buoni match cmq! (?)
  M[M<0] = 0
  A = matrix(rep(NA,length(M)), nrow = nrow(M) )
  best=list(row=rep(0,nrow(M)),col=rep(0,nrow(M)),similarity=rep(0,nrow(M)))
  for (i in 1:nrow(M)){
    for(j in 1:ncol(M)){
      A[i,j] = max( max(A[i-1,j],0), max(A[i,j-1],0), max(A[i-1,j-1],0)+M[i,j] ,na.rm = T)
    }
  }
  A2 = A
  while(sum(A2)){
    x = which(A2 == max(A2), arr.ind = TRUE)[1,]
    A2[x[1]:nrow(A2),] = 0
    A2[,x[2]:ncol(A2)] = 0
    best$row = c(best$row,x[1])
    best$col = c(best$col,x[2])
    best$similarity = c(best$similarity, M[x[1],x[2]])
  }
  xbest = data.frame(best)
  xbest = xbest[order(xbest$row),]
  xbest = xbest[xbest$col!=0 & xbest$row !=0,]
  xbest$lag =   allpik2[xbest$col]-allpik1[xbest$row]
  xbest$s1 = allpik1[xbest$row]
  xbest$s2 = allpik2[xbest$col]
  #instantiate new sync class object
  signal$PMBest = PMBest(NULL,NULL,xbest,lagSec,sgol_p,sgol_n,weightMalus)
  signal
}


ppSync = function(signal, type=c("block","continuous")) {
  if(is.null(signal$PMBest)||is.null(signal$PMBest$xBest))stop("Please run ppBest before.")
  type=match.arg(type)
  d = signal$s2
  d2  = signal$s1
  lagSec = attr(signal$PMBest, "lagSec")
  sampRate = sampRate(signal)
  ransamp = lagSec * sampRate
  if(is.null(signal$ccf$ppBest)) stop("you need to run ppBest beforehand")
  xbest = signal$PMBest$xBest
  xbest$sync = rep(NA,nrow(xbest))
  
  lagvec = rep(NA,length(d) )
  lags = xbest$lag
  for(i in 1:nrow(xbest)){
    #per il primo picco porta tutti i lag al valore del primo picco
    if(i==1) {lagvec[1:xbest$s1[1]] = lags[1]
    } else {
      # per tutti gli altri lag fai una interpolazione
      ab = (xbest$s1[i-1]+1):xbest$s1[i]
      lagvec[ab] = round(seq(from=lags[i-1], to=lags[i], length.out = length(ab)))
    }
  }
  lagvec[is.na(lagvec)] = lags[length(lags)]#fix the last values
  
  if(type=="block"){
    syncvec = rep(NA,length(d))
    for(i in 1:(nrow(xbest)-1)){
      ## trova il segmento più corto e allungalo fino alla larghezza di quello di più lungo
      ## mantenendolo al centro
      ab = list()
      ab[[1]] = xbest$s1[i]:xbest$s1[i+1]
      ab[[2]] = xbest$s2[i]:xbest$s2[i+1]
      short = s = which.min(lapply(ab,length))
      long =  l = which.max(lapply(ab,length))
      
      ### [MODE 1] use the width of the longer, by enlarging the WINDOW of the shorter #########
      delta = abs(length(ab[[1]]) - length(ab[[2]]))
      if(delta%%2!=0){ #se delta è dispari, invece di usare delta/2 ti tocca fare delta-1 /2 e rimettere l'1 alla fine
        ab[[short]] = (min(ab[[short]]) - (delta-1)/2) : (max(ab[[short]]+ (delta-1)/2) +1)
      } else { #se pari o uguale a zero (stessa lunghezza)
        ab[[short]] = (min(ab[[short]]) - delta/2) : (max(ab[[short]]+ delta/2))
      }
      ## controllo che non si vada oltre ai boundaries 1:length(d)
      toKeep1 = which(ab[[1]]>0 & ab[[1]]<length(d))
      toKeep2 = which(ab[[2]]>0 & ab[[2]]<length(d2))
      ab = lapply(ab, function(x) x[if(length(toKeep1)<length(toKeep2)) toKeep1 else toKeep2  ])
      iCor = cor(d[ab[[1]]] , d2[ab[[2]]]) #"pear"
      
      ### [MODE 2] stretch the shorter to the width of the longer #########
      # #is this working?
      # data = list("s1" = d, "s2" = d2)
      # 
      # if(length(ab[[1]]) != length(ab[[2]])){
      #   toCorL = data[[l]][ab[[l]]]
      #   toCorS = approx(x=1:length(ab[[s]]),y=data[[s]][ab[[s]]],
      #                   xout=seq(1,length(ab[[s]]),length.out = length(ab[[l]])))$y
      # } else {
      #   toCorL = data[[1]][ab[[1]]]
      #   toCorS = data[[2]][ab[[2]]]
      # }
      # iCor = cor(toCorL,toCorS)
      
      ### [FINALLY] save objects ##########
      xbest$sync[i] = iCor
      syncvec[xbest$s1[i]:xbest$s1[i+1]] = iCor
    }
    
  } else if (type == "continuous") {
    ############################################################
    ## now generate the interpolated bestlag and bestccf
    stop("not yet developed")
    # winSec = signal$ccf$settings$winSec
    # incSec = signal$ccf$settings$incSec
    # win = winSec*sampRate
    # lagSamp = lagSec*sampRate
    # inc = incSec * sampRate
    # n_win = ceiling((length(d)-win-lagSamp+1)/inc) #calcola il numero di finestre per ciascun file
    # 
    # ## adesso trova i valori di lag corrispondenti per ciascun inc, e calcola la correlazione
    # l_lag = lagvec[win/2 + (seq_len(n_win)-1)*inc]  # valore di lag nella posizione centrale di ciascuna finestra
    # a_cli = 1 + lagSamp + (seq_len(n_win)-1)*inc    # sample iniziale finestra del Clinician
    # b_cli = a_cli + win -1                              # sample finale finestra del Clinician
    # a_pat = a_cli - l_lag                           # sample iniziale finestra del Patient
    # b_pat = b_cli - l_lag                           # sample finale finestra del Patient
    #                                                 # Valori positivi di lag -> patient leading
    # 
    # l_cor=sapply(seq_len(n_win), function(iWin) { #-----per ciascuna finestra--------
    #   x = d [a_cli[iWin]:b_cli[iWin]]                     # estrai finestra del Clinician
    #   y = d2[a_pat[iWin]:b_pat[iWin]]             # e del patient. Valori positivi di lag -> patient leading.
    #   if(all(is.na(x)) || all(is.na(y))) NA
    #   else  suppressWarnings(cor(x,y,method="pearson",use="complete.obs"))
    # })
    # 
    # res = data.frame( "ppCor"= l_cor, "ppLag"=round(l_lag/sampRate,2), "pp_pat_sam" = a_pat + win/2, "pp_clin_sam" = a_cli + win/2)
    # 
    # # if(length(signal$ccf$ccfmat)==0) {
    # #   signal$ccf$ccfmat = res
    # # } else signal$ccf$ccfmat = cbind(signal$ccf$ccfmat, res)
  }
  
  #qui salva l'output
  signal$PMBest$xBest = xbest
  signal$PMBest$sync = DyadStream(syncvec, "PMBest_Sync", sampRate, start=start(d), col=rgb(191,50,59,max=255))
  signal$PMBest$lag = DyadStream(lagvec, "PMBest_Lag", sampRate, start=start(d), col=rgb(253,177,2,max=255))
  return(signal)
  
}

##peakfinder è una funzione ausiliaria di ppBestLag, è sviluppata in "peak-centered-flex-CCF_v1.93"
peakFinder = function(x, sgol_p = 2, sgol_n = 25, mode=c("peaks","valleys"), correctionRangeSeconds = 0.5, valid){
  if(missing(valid)) valid = rep(T,length(x))
  sampRate = frequency(x)
  fdderiv1  = ts(signal::sgolayfilt(x,  p =sgol_p, n = sgol_n, m=1),frequency = 10,start=start(x))
  mode = match.arg(mode)
  pik = sign(embed(fdderiv1,2)) #embed appaia al segnale il segnale laggato di 1 samp
  if(mode=="peaks")
    pikboo = c(pik[,2] - pik[,1] ==  2, FALSE) #embed perde 1 sample, quindi aggiungi un FALSE alla fine
  else
    pikboo = c(pik[,2] - pik[,1] == -2, FALSE) #embed perde 1 sample, quindi aggiungi un FALSE alla fine
  
  pikboo[valid==F] = FALSE #cancella i picchi nelle aree con artefatti
  
  piksam = which(pikboo) #in quali sample c'è un'inversione di segno della derivata?
  #correzione manuale: cerca il valore più alto nei dintorni del cambio di derivata
  for(v in seq_along(piksam)){
    #individua il range con 0.5s prima e dopo la valle della derivata (esclusi gli estremi inzio e fine ts)
    search_interval = x[max(1,piksam[v]-correctionRangeSeconds*sampRate):min(length(x),piksam[v]+correctionRangeSeconds*sampRate)]
    #trova il più piccolo e aggiorna pikmsam
    if(mode=="peaks")
      piksam[v] = max(1, piksam[v]+  (which.max(search_interval) - round(length(search_interval)/2)))
    else
      piksam[v] = max(1, piksam[v]+  (which.min(search_interval) - round(length(search_interval)/2)))
    piksam[v] = min(piksam[v], length(x))
  }
  #ricrea pikmboo dai sample corretti
  pikboo = rep(F,length(pikboo))
  pikboo[piksam] = T
  piks = time(x)[piksam]
  list("bool" = pikboo,
       "samples" = piksam,
       "seconds" = piks)
}

