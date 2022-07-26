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




#' pmBest
#'
#' @param experiment 
#' @param signals 
#' @param lagSec 
#' @param sgol_p 
#' @param sgol_n 
#' @param weightMalus weightmalus è la percentuale di malus per il lag più estremo. Es, con weightMalus= 20 e r = 1 al massimo lag, la trasformazione diventa r' = 0.8
#' @param match_threshold a value between 0 and 1, specifying the similarity threshold to assign a peak-peak match.
#' @param minSizeSec the smallest allowed window size. very small windows may artificially inflate synchrony.
#' @param algorithm
#' @param outputName
#' @param scaledCorrelation logical. If TRUE, correlations will be scaled using the whole signals range, and be much more robust
#' @param correctionRangeSeconds parameter passed to peakFinder, seconds of manual correction for derivative based peak-finding
#' @param minPeakDelta parameter passed to peakFinder, minimum valley-peak delta to detect a peak
#' @return
#' @export
#'
#' @examples
pmBest = function(experiment, signals="all", lagSec=20,
                  sgol_p = 2, sgol_n = 25,  weightMalus = 20,
                  match_threshold = 0.25,minSizeSec=1, algorithm=c("classic","AMICo","sccf"), outputName = "PMBest",scaledCorrelation=F, 
                  correctionRangeSeconds = 0.5, minPeakDelta){
  algorithm=match.arg(algorithm)
  if(!is(experiment,"DyadExperiment")) stop("Only objects of class DyadExperiment can be processed by this function")
  cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]),
             " follows the ", s1Name(experiment[[1]]),"\r\n")) #verified!
  nSessions = length(experiment)
  experiment2 = Map(function(session,iSession){
    if(signals=="all") signals = names(session)
    cat("\r\n",paste(dyadId(session),sessionId(session)))
    session[signals] = Map(function(signal){
      cat(" |",name(signal))
      #the first function calculates best lag
      signal = peakMatch(signal, lagSec=lagSec, sgol_p=sgol_p, sgol_n=sgol_n,
                         weightMalus=weightMalus, match_threshold=match_threshold, outputName=outputName,correctionRangeSeconds=correctionRangeSeconds, minPeakDelta=minPeakDelta)
      #the second calculates the actual sync values
      if(algorithm == "AMICo")
        signal = ppSync_dev(signal,minSizeSec,outputName=outputName, scaledCorrelation)
      else stop("other algorithms are deprecated")
      # else if(algorithm == "classic")
      #   signal = ppSync(signal,minSizeSec,outputName=outputName)
      # else if(algorithm == "sccf")
      #   signal = ppSync_sccf(signal,outputName=outputName)
      
      return(signal)
    }, session[signals])
    return(session)
  },experiment,seq_along(experiment))
  cat("\r\nDone ;)\r\n")
  classAttr(experiment2)=classAttr(experiment)
  return(experiment2)
}


## peak picking best lag
peakMatch = function(signal,lagSec=4, sgol_p = 2, sgol_n = 25, weightMalus = 30,match_threshold=0.0, outputName, correctionRangeSeconds, minPeakDelta) {
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  
  #signal = lr$all_vid1_01$contempt
  # signal = lr10$all_CC_1$SC
  # signal = lr10f$all_02_1$SC
  #pacs best settings
  ## lr = pmBest(d,signals = "all",lagSec=4,match_threshold=0.25,
  #             minSizeSec=5,weightMalus = 35 ,algorithm = "AMICo",outputName = "PMdev")
  
  d = signal$s1
  d2  = signal$s2
  sampRate = frequency(signal)
  ransamp = lagSec * sampRate
  
  
  ### peaks-valleys detection ########
  s1b = peakFinder(d,  sgol_p, sgol_n, mode = "b", correctionRangeSeconds, minPeakDelta) 
  s2b = peakFinder(d2, sgol_p, sgol_n, mode = "b", correctionRangeSeconds, minPeakDelta)

  allSec1 = time(d)[s1b$samples]
  allSec2 = time(d2)[s2b$samples]
  
  #full matrix with all matches values
  M = matrix(NA ,nrow=length(s1b$samples),ncol=length(s2b$samples))
  
  for(i in 1:length(s1b$samples)){
    ipik = s1b$samples[i] 
    #trova il range attorno al picco i in cui cercare la lag
    ###NB nella v1.1 questo avviene da valle a valle, per la versione fissa guarda v1.0c
    search_range = (ipik-ransamp):(ipik+ransamp)
    
    ## seleziona il range di confronto va dal picco/valle precedente a quello successivo di ipik
    #se non ci sono picchi/valli prima, parti dall'inizio del segnale
    a = ifelse(any(s1b$samples<ipik),  max(s1b$samples[s1b$samples<ipik]), 1)
    #se non ci sono picchi/valli dopo, usa la fine del segnale
    b = ifelse(any(s1b$samples>ipik),  min(s1b$samples[s1b$samples>ipik]), length(d))
    ab = a:b
    
    search_range[search_range<=0] = NA #this should not be needed...
    ab[ab<=0] = NA #this should not be needed...
    matches = which(s2b$bool[search_range]) # trova picchi in d2 nell'intorno di ipik (su d)
    matches = matches + (ipik-ransamp) -1     # passa da posizione relativa a search_range a posizione assoluta su d2
    matches_i = which(s2b$samples %in% matches)   
    nMatch = length(matches)
    
    if(nMatch>0) {
      ###per ogni match
      for(f in 1:nMatch){
        ipik2 = matches[f]
        lagv = matches[f] - ipik #distanza fra i due picchi, in samples. Valori negativi indicano giallo anticipa blu
        
        #trova il range valle-valle o picco-picco del segno matchato
        a2 = ifelse(any(s2b$samples<ipik2),  max(s2b$samples[s2b$samples<ipik2]), 1)
        b2 = ifelse(any(s2b$samples>ipik2),  min(s2b$samples[s2b$samples>ipik2]), length(d2))
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
        
        toCor1 = d[range1]
        toCor2 = d2[range2]
        dd1 = as.numeric(diff(toCor1)); dd2 = as.numeric(diff(toCor2))
        
        # #invece delle correlazioni usa la differenza fra le derivate normalizzate
        #questo è la formula classica di Agosto 2018 e precedenti. Gold standard
        m2d = 1-mean(abs(rangeRescale(dd1,0,1) - rangeRescale(dd2,0,1)))
        
        # #per favore normalizzare le derivate così non ha senso! mantieni almeno il segno:
        # dd1x = rangeRescale(dd1,-1,1,pres.signs = T)
        # dd2x = rangeRescale(dd2,-1,1,pres.signs = T)
        # distance = sqrt(abs(dd1x - dd2x)) # radice quadrata serve a normalizzare la distribuzione
        # m2d = sqrt(2) - distance
        
        thisCorOld = mean(m2d)^2
        # thisCorOld = mean(m2d^2)
        
        
        
        # marciCor = rangeRescale(cor(diff(toCor1),diff(toCor2)),0,1,-1,1)
        # 
        # ## se normalizzi ciascuna derivata, la differenza si riduce arbitrariamente.
        # ##prova a normalizzare sulla base di minimi e massimi di entrambe le serie:
        # 
        # dmin = min(dd1,dd2,na.rm = T); dmax = max(dd1,dd2,na.rm = T)
        # dd1 = rangeRescale(dd1,0,1,dmin,dmax)
        # dd2 = rangeRescale(dd2,0,1,dmin,dmax)
        # 
        # #•calcola l'area di differenza massima teorica:
        # n=length(dd1)
        # yy = abs(sin(seq(0,2*pi,length.out = n))/2)
        # amax = (2*sum(yy))
        # m2d = amax-sum(abs(dd1-dd2))
        # # m2d = 1-sqrt(mean((dd1-dd2)^2))# m2d = 1-mean(abs(dd1-dd2))
        # thisCor = m2d/n
        # ## bel tentativo ma farlocco. Infatti la distanza media fra le derivate è sempre più
        # ## o meno uguale. è più interessante la correlazione, o la distanza delle derivate
        # ## normalizzate.
        # 
        # #ultimo tentativo! Normalizza le slope all'inizio, quindi 0.5 è il massimo individuale!
        # ## poi fai 1- la distanza assoluta fra le slope (max=1), che da la similitudine
        # ## tuttavia, se il segno è diverso, qualsiasi valore di "similitudine" in realtà
        # ## indica un trend diverso nel segnale, quindi assegna segno meno
        # dd1 = sd1[range1]
        # dd2 = sd2[range2]
        # distanza = abs(dd1-dd2)
        # simil = 1 - abs(dd1-dd2)
        # res = mean((sign(dd1)*sign(dd2))*simil)
        # res = rangeRescale(res,0,1,-1,1)
        # ## Non male, ma thisCorOld resta il migliore lol!
        
        thisCor = thisCorOld
        
        
        ## plot area ##########################
        # tit = paste0("F:",f," - dd:",round(thisCor,2)," | sCor:",round(marciCor,2)," | ddOld:",round(thisCorOld,2)," | res:",round(res,2))
        # 
        # plot(range1[-1],yy+0.5,ylim=c(0,1),type="h",main=tit,col="grey80");lines(range1[-1],-yy+0.5,type="h",col=0)
        # abline(h=0.5)
        # lines(ab,        rangeRescale(d[ab],  0,1,xmin = min(d[ab],d2[ab2]),xmax = max(d[ab],d2[ab2])),col="blue");
        # lines(ab2 - lagv,rangeRescale(d2[ab2],0,1,xmin = min(d[ab],d2[ab2]),xmax = max(d[ab],d2[ab2])),col="red")
        # 
        # lines(range1,dd1+0.5,col="grey30",lwd=2);
        # lines(range2-lagv,dd2+0.5,col="grey1",lwd=2);
        
        # ###################################
        
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
  ## ignora le correlazioni sotto un certo threshold:
  M[M<match_threshold] = 0
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
  xbest$lag =   s2b$samples[xbest$col]-s1b$samples[xbest$row]
  xbest$s1 = s1b$samples[xbest$row]
  xbest$s2 = s2b$samples[xbest$col]
  #instantiate new sync class object
  signal[[outputName]] = PMBest(NULL,NULL,xbest,lagSec,sgol_p,sgol_n,weightMalus)
  signal
}


ppSync = function(signal,minSizeSec, outputName) {
  if(is.null(signal[[outputName]])||is.null(signal[[outputName]]$xBest))stop("Please run peakMatch before.")
  d = signal$s1
  d2  = signal$s2
  lagSec = attr(signal[[outputName]], "lagSec")
  sampRate = sampRate(signal)
  ransamp = lagSec * sampRate
  xbest = signal[[outputName]]$xBest
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
  lagvec[is.na(lagvec)] = lags[length(lags)] #fix the last values
  
  
  syncvec = rep(NA,length(d))
  abCum = list()
  iCum = c()
  deltaVec = rep(NA,nrow(xbest))
  xbest$syncBound = c(T,rep(F,nrow(xbest)-1))
  xbest$syncStart = rep(F,nrow(xbest))
  xbest$syncEnd   = rep(F,nrow(xbest))
  
  for(i in 1:(nrow(xbest)-1)){
    
    ab = list()
    ab[[1]] = xbest$s1[i]:xbest$s1[i+1]
    ab[[2]] = xbest$s2[i]:xbest$s2[i+1]
    ii = i
    
    if(length(abCum)>0){ #se ci sono dei dati lasciati indietro aggiungili ad ab
      ab[[1]] = c(abCum[[1]],ab[[1]])
      ab[[2]] = c(abCum[[2]],ab[[2]])
      abCum = list() # e resetta abCum
      ii = c(iCum, ii)
      iCum = c()
    }
    #if the longest of ab is at least minSizeSec length
    if(length(ab[[which.max(lapply(ab,length))]]) >= minSizeSec*sampRate){
      ## trova il segmento più corto e allungalo fino alla larghezza di quello di più lungo
      ## mantenendolo al centro
      short = s = which.min(lapply(ab,length))
      long =  l = which.max(lapply(ab,length))
      
      ### [MODE 1] use the width of the longer, by enlarging the WINDOW of the shorter #########
      delta = abs(length(ab[[1]]) - length(ab[[2]]))
      deltaVec[i]=delta
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
      # iCor = cor(diff(d[ab[[1]]]) , diff(d2[ab[[2]]])) #"pear" on slopes
      
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
      
      
      ### diagnostic plot ####
      # rs1 = rangeRescale(d,-1,1)
      # rs2 = rangeRescale(d2,-1,1)
      # par(mfrow=c(1,2))
      # plot(rs1[ab[[1]]],type="l",main=iCor,ylim=c(-1,1))
      # lines(rs2[ab[[2]]])
      # #scatter
      # plot(d[ab[[1]]],d2[ab[[2]]]);abline(lm(d[ab[[2]]]~d2[ab[[1]]]),lty=3);
      # midLine(mean(d[ab[[1]]]),mean(d2[ab[[2]]]))
      # midLine = function(ux,uy){lines(seq(-1+ux,1+ux,length.out = 100),seq(-1+uy,1+uy,length.out = 100))}
      
      
      ### [FINALLY] save objects ##########
      xbest$sync[ii] = iCor
      xbest$syncBound[i+1] = T
      #start of each sync windows is located in the center between begin of s1 and s2 windows
      xStart = round(mean(min(ab[[1]]),min(ab[[2]])))
      xEnd = round(mean(max(ab[[1]]),max(ab[[2]])))
      syncvec[xStart:xEnd] = iCor
      xbest$syncStart[i] = xStart
      xbest$syncEnd[i] = xEnd
      
    }  else { #salva i dati
      abCum = ab
      iCum = ii
    }
  } # end of for
  
  #qui salva l'output ############
  syncvec = ts(syncvec, start=start(d), frequency = sampRate)
  lagvec  = ts(lagvec,  start=start(d), frequency = sampRate)
   # applica gli artefatti
  for(i in seq_len(nrow(signal$artefacts))){
    window(syncvec,signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
    window(lagvec,signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
  }
  signal[[outputName]]$xBest = xbest
  signal[[outputName]]$sync = DyadStream(syncvec, "PMBest_Sync", col=rgb(191,50,59,max=255))
  signal[[outputName]]$lag = DyadStream(lagvec, "PMBest_Lag", col=rgb(253,177,2,max=255))
  return(signal)
  
}

#AMICo: Adaptive Matching Interpolated Correlation
ppSync_dev = function(signal,minSizeSec, outputName, scaledCorrelation) {
  ##questa è una versione di debug e sviluppo di ppSync di rIP
  # richiede dei dati generati con rIP in un oggetto chiamato lr
  # ottimi quelli di "PACS_rIP engine.R"
  #
  #la versione attuale (l'ultima su github ) è imperfetta:
  # 1. usare la finestra più lunga causa degli errori in caso di cambio di lag:
  #   ^  ^
  # /  \/ \
  #   ^ ^    <-qui anche se i due picchi sono uguali la cor viene bassa
  # meglio se ALMENO una finestra sia sufficientemente lunga e quella breve venga interpolata
  #
  # 2. aggregando finestre diverse per avere una lunghezza minima si sminchia l'appaiamento fra picchi
  # cosa che si vede già con ii = c(2,3) {II e III righe di xbest} che vengono aggregate ma il lag risultante è errato.
  
  ## idee: - pesare per la differenza di durata?
  ##       - pesare per la differenza di ampiezza (normalizzata sulla SD individuale?)
  
  cat(" - AMICo algorithm v.1.1")
  #please refer to ppSync dev.R in research folder of DyadClass
  
  # signal = lr$all_CC_1$SC
  # minSizeSec=5 <-- can be smaller as the check is on the shorter segment?
  if(is.null(signal[[outputName]])||is.null(signal[[outputName]]$xBest))stop("Please run peakMatch before.")
  if(nrow(signal[[outputName]]$xBest)<2) stop("peakMatch found only one peak in the signal")
  if(nrow(signal[[outputName]]$xBest)<10) warning("peakMatch found less than 10 peaks")
  d = signal$s1
  d2  = signal$s2
  lagSec = attr(signal[[outputName]], "lagSec")
  sampRate = sampRate(signal)
  ransamp = lagSec * sampRate
  xbest = signal[[outputName]]$xBest
  xbest$sync = rep(NA,nrow(xbest))
  
  #crea il vettore di lag al sample rate finale
  lagvec = rep(NA,length(d) )
  lags = xbest$lag
  for(i in 1:nrow(xbest)){
    #per il primo picco porta tutti i lag al valore del primo picco
    if(i==1) {lagvec[1:xbest$s1[1]] = lags[1]
    } else {
      # per tutti gli altri lag fai una interpolazione
      #così da avere dei valori progressivi di lag da xbest i-1 a i
      ab = (xbest$s1[i-1]+1):xbest$s1[i]
      lagvec[ab] = round(seq(from=lags[i-1], to=lags[i], length.out = length(ab)))
    }
  }
  lagvec[is.na(lagvec)] = lags[length(lags)] #fix the last values
  
  #crea il vettore di sync al sample rate finale
  syncvec = rep(NA,length(d))
  iCum = c()
  deltaVec = rep(NA,nrow(xbest))
  xbest$syncBound = c(T,rep(F,nrow(xbest)-1))
  xbest$syncStart = rep(F,nrow(xbest))
  xbest$syncEnd   = rep(F,nrow(xbest))
  data = list("s1" = d, "s2" = d2)
  # datax = list("s1" = rangeRescale(d/sd(d),0,1), "s2" = rangeRescale(d2/sd(d2),0,1)) #vedi cosa cambia normalizzando il segnale per l'SD di tutta la seduta
  # 
  #for each match (i.e. xbest row)
  
  #alpha normalizing elements:
  # nMin = min(quantile(d,probs = 0.05), quantile(d2,probs=0.05)) #smokin
  # nMax = max(quantile(d,probs = 0.95), quantile(d2,probs=0.95)) #shit
  
  dd2 = c(d,d2)
  normPam = sd(dd2)
  dd2n = scale(dd2)
  maxDist = sqrt((max(c(d,d2)) - min(c(d,d2)))^2)
  
  for(i in 1:(nrow(xbest)-1)){
    # i=0
    # i=i+1
    iCum = c(iCum,i)
    # ii = i
    #define working samples interval
    ab = list()
    ab[[1]] = xbest$s1[i]:xbest$s1[i+1]
    ab[[2]] = xbest$s2[i]:xbest$s2[i+1]
    short = s = which.min(lapply(ab,length))
    long =  l = if(s==2) 1 else 2
    
    toCorL = data[[l]][ab[[l]]]
    if(any(is.na(data[[s]][ab[[s]]])) || is.null(data[[l]][ab[[l]]]) || length(data[[l]][ab[[l]]])==0){
      cat("\r\nWarning: NAs found!\r\nFile: ",name(signal),"\r\nxbest row: ",i,"\r\n" ,data[[s]][ab[[s]]])
    }
    if(length(ab[[s]])<2) {
      ab[[s]] = rep(ab[[s]],2)
      cat(" - Warning xbest repeated peak at line ", i)
    } #this is a bad hack for when xbest connects twice to the same peak.
    toCorS = approx(x=1:length(ab[[s]]),y=data[[s]][ab[[s]]],                    ##this works fine!
                    xout=seq(1,length(ab[[s]]),length.out = length(ab[[l]])))$y  ##
    # toCorLx = datax[[l]][ab[[l]]]
    # toCorSx = approx(x=1:length(ab[[s]]),y=datax[[s]][ab[[s]]],                    ##this works fine!
    #                  xout=seq(1,length(ab[[s]]),length.out = length(ab[[l]])))$y  ##
    
    
    rs1 = if(l==1) toCorL else toCorS
    rs2 = if(l==2) toCorL else toCorS
    # rs1x = if(l==1) toCorLx else toCorSx
    # rs2x = if(l==2) toCorLx else toCorSx
    
    #if previous intervals were to small abCum will be > 0
    if(length(iCum)>1){ #se ci sono dei dati lasciati indietro anteponili ad rs1/rs2
      rs1 = c(rs1mem,rs1)
      rs2 = c(rs2mem,rs2)
      # rs1x = c(rs1memx,rs1x)
      # rs2x = c(rs2memx,rs2x)
    }
    
    if(length(rs1) > minSizeSec*sampRate){
      # #normalizing element
      if(scaledCorrelation){
        # #his is shit
        # rs1=c(nMin,nMax,rs1)
        # rs2=c(nMin,nMax,rs2)
        rs1 = rs1-mean(rs1)
        rs1 = rs1/normPam
        rs2 = rs2-mean(rs2)
        rs2 = rs2/normPam
        
        iCor = sqrt(sum((rs2-rs1)^2))
        iCor = 1- iCor/maxDist
      }else 
      iCor = cor(rs1,rs2, use = "c")
      
      # ### diagnostic plot ####
      # par(mfrow=c(2,2))
      # 
      # # what the correlation sees:
      # rs1p = rangeRescale(rs1,0,1)
      # rs2p = rangeRescale(rs2,0,1)
      # plot(rs1p,type="l",main=paste("correlation sees | r:",round(iCor,2)),ylim=c(0,1), col=2)
      # lines(rs2p, col=4)
      # ##what I see:
      # s1range = xbest$s1[iCum[1]]:xbest$s1[iCum[length(iCum)]+1]
      # s2range = xbest$s2[iCum[1]]:xbest$s2[iCum[length(iCum)]+1]
      # s1y= d[s1range]; s2y = d2[s2range]
      # 
      # plot(s1range,s1y,
      #      xlim=c(min(s1range,s2range),max(s1range,s2range)), ylim=c(min(s1y,s2y),max(s1y,s2y)),
      #      type="l",main=paste("i see | r:",round(iCor,2)), col=2)
      # lines(s2range,s2y, col=4)
      # 
      # #normalizing element
      # rs1=c(nMin,nMax,rs1)
      # rs2=c(nMin,nMax,rs2)
      # iCor = cor(rs1,rs2, use = "c")
      # # iCor = cor(diff(rs1),diff(rs2), use = "c")
      # 
      # # what the correlation sees:
      # rs1p = rangeRescale(rs1,0,1)
      # rs2p = rangeRescale(rs2,0,1)
      # plot(rs1p,type="l",main=paste("correlation sees | r:",round(iCor,2)),ylim=c(0,1), col=2)
      # lines(rs2p, col=4)
      # ##what I see:
      # s1range = xbest$s1[iCum[1]]:xbest$s1[iCum[length(iCum)]+1]
      # s2range = xbest$s2[iCum[1]]:xbest$s2[iCum[length(iCum)]+1]
      # s1y= d[s1range]; s2y = d2[s2range]
      # s1y=c(nMin,nMax,s1y)
      # s2y=c(nMin,nMax,s2y)
      # s1range = c(min(s1range)-1,min(s1range)-2,s1range)
      # s2range = c(min(s2range)-1,min(s2range)-2,s2range)
      # plot(s1range,s1y,
      #      xlim=c(min(s1range,s2range),max(s1range,s2range)), ylim=c(min(s1y,s2y),max(s1y,s2y)),
      #      type="l",main=paste("i see | r:",round(iCor,2)), col=2)
      # lines(s2range,s2y, col=4)

      # ##dividing by whole signal sd
      # plot(rs1x,type="l",main=cor(rs1x,rs2x),ylim=c(0,1), col=2)
      # lines(rs2x, col=4)
      # ## instead of correlation you may know calculate a dissimilarity index based on standard deviations differences.
      # ## eg sum(abs(rs1x - rs2x))/n
      # ## but it has problems:
      # ## 1. SD on nonstationary data is bad. should you use local SD? does it still make sense?
      # ## 2. how do you normalize it? there is no theoretical maximum! sd(s1)*sd(s2) what does it means?
      # 
      # #scatter
      # plot(rs1-rs2);abline(lm(d[ab[[2]]]~d2[ab[[1]]]),lty=3);
      # midLine(mean(d[ab[[1]]]),mean(d2[ab[[2]]]))
      # midLine = function(ux,uy){lines(seq(-1+ux,1+ux,length.out = 100),seq(-1+uy,1+uy,length.out = 100))}
      # 
      # 
      
      ### [FINALLY] save objects ##########
      xbest$sync[iCum] = iCor
      xbest$syncBound[i+1] = T
      #start of each sync windows is located in the center between begin of s1 and s2 windows
      xStart = round(mean(c(xbest$s1[iCum[1]],xbest$s2[iCum[1]] )))+1
      xEnd = round(mean(c(xbest$s1[iCum[length(iCum)]+1],xbest$s2[iCum[length(iCum)]+1] )))
      all(is.na(syncvec[xStart:xEnd]))
      syncvec[xStart:xEnd] = iCor
      xbest$syncStart[i] = xStart
      xbest$syncEnd[i] = xEnd
      
      iCum = NULL #reset iCum
      
    } else {
      rs1mem = rs1
      rs2mem = rs2
      # rs1memx = rs1x
      # rs2memx = rs2x
    }
  } # end of for
  
  #qui salva l'output ############
  syncvec = ts(syncvec, start=start(d), frequency = sampRate)
  lagvec  = ts(lagvec,  start=start(d), frequency = sampRate)
  # applica gli artefatti
  for(i in seq_len(nrow(signal$artefacts)) ){ 
    #@TSBUG
    #in alcune installazioni di R c'è un bug per cui window(x xstart(x), xend(x)) risulta 1 sample più lungo di x
    if(signal$artefacts[i,"end"] != signal$artefacts[i,"start"]){
      
      realEnd = c(signal$artefacts[i,"end"]-1, frequency(syncvec))
      realStart = c(signal$artefacts[i,"start"], 1)
      window(syncvec,realStart, realEnd ) <- NA
      window(lagvec,realStart, realEnd ) <- NA
    }
  }
  signal[[outputName]]$xBest = xbest
  signal[[outputName]]$sync = DyadStream(syncvec, "PMBest_Sync", col=rgb(191,50,59,max=255))
  signal[[outputName]]$lag = DyadStream(lagvec, "PMBest_Lag", col=rgb(253,177,2,max=255))
  return(signal)
}

ppSync_sccf = function(signal, winSec = 10, incSec=1, outputName) {
  ##questa è una versione di debug e sviluppo di ppSync di rIP
  # le altre due versioni di ppSync hanno un approccio molto deterministico
  # calcolando la correlazione esattamente di ciascun picco con ciascun picco
  # tra le altre cose ha lo svantaggio di non avere un valore continuo
  # inoltre la performance è misera confrontata con lag0
  # qua invece uso la lag calcolata in xbest per muovere le finestre di cross-correlazione
  # risultando una specie di ibrido fra CCFBest e PMBest
  # ottenendo però valori continui basati su lag precisa.
  
  #import C correlation function
  C_cor=get("C_cor", asNamespace("stats"))
  cat(" - wccs version of ppSync")
  #please refer to ppSync dev.R in research folder of DyadClass
  
  # signal = mimic$s09_01$SC
  # outputName = "PMsccf"
  # minSizeSec=5 <-- can be smaller as the check is on the shorter segment?
  if(is.null(signal[[outputName]])||is.null(signal[[outputName]]$xBest)) stop("Please run peakMatch before.")
  
  lagSec = attr(signal[[outputName]], "lagSec")
  sampRate = sampRate(signal)
  ransamp = lagSec * sampRate
  xbest = signal[[outputName]]$xBest
  # xbest$sync = rep(NA,nrow(xbest))
  
  iFile = data.frame(movAv(diff(signal$s1),5,1),
                     movAv(diff(signal$s2),5,1))
  
  #crea il vettore di lag al sample rate finale
  lagvec = rep(NA,nrow(iFile) )
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
  lagvec[is.na(lagvec)] = lags[length(lags)] #fix the last values
  # lagvecsec = lagvec[seq(1,length(lagvec),by=sampRate)]
  
  #ora CCF!!
  win = winSec*sampRate
  lagSamp = lagSec*sampRate
  inc = incSec * sampRate
  
  n_win = ceiling((nrow(iFile)-win-2*lagSamp+1)/inc) #calcola il numero di finestre
  lcc=sapply(seq_len(n_win)-1, function(iWin)
  { #-----per ciascuna finestra--------
    # ab è il range di sample di ciascuna finestra, partendo da lagSamp invece che da 1
    # in modo da non avere numeri negativi applicando il lag
    ab = ((iWin*inc +1):(iWin*inc +win))+lagSamp 
    iLag = lagvec[(iWin*inc +1)+win/2+lagSamp] #valore di lag al centro della finestra
    # plot(ab,rangeRescale(signal$s1[ab],1,0),type="l");lines(ab,rangeRescale(signal$s2[ab],1,0),lty=3);lines(ab,rangeRescale(signal$s2[ab+iLag],1,0),lty=3,col="red")
    xWin = iFile[ab,1];
    yWin = iFile[ab+iLag,2] #estrai i dati della finestra in un vettore per sogg x e y
    # plot(ab,rangeRescale(iFile[ab,1],1,0),type="l");lines(ab,rangeRescale(iFile[ab+iLag,2],1,0),lty=3)
    if(sum(xWin!=yWin,na.rm = T)<length(xWin)/2 )  NA #controlla che ci siano almeno 3 punti !=0
    else {
      suppressWarnings(.Call(C_cor, xWin, yWin, 2L, FALSE))
    }
    
    
  }) #fine lapply finestre
  
  
  #qui salva l'output ############
  signal[[outputName]]$xBest = xbest
  xStart = c(start(signal)[1] +round(winSec/2)+lagSec,start(signal)[2])

  signal[[outputName]]$sync = DyadStream(lcc, "PMBest_Sync", col=rgb(191,50,59,max=255), start=xStart, frequency=1/incSec )
  signal[[outputName]]$lag = DyadStream(lagvec, "PMBest_Lag", col=rgb(253,177,2,max=255), start=start(signal),frequency=sampRate )
  
  # applica gli artefatti
  warning("artefact removal in ppSync_sccf is untested. window.DyadStream '<-' behaviour is unknown")
  for(i in seq_len(nrow(signal$artefacts))){
    window(signal[[outputName]]$sync,signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
    window(signal[[outputName]]$lag,signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
  }
  return(signal)
  
}

##peakfinder è una funzione ausiliaria, è sviluppata in "peak-centered-flex-CCF_v1.93"
#' Title
#'
#' @param x a ts object or a numeric vector.
#' @param sgol_p smoothing filter order. If NA, the filtering is disabled.
#' @param sgol_n smoothing filter length (must be odd).
#' @param mode should the function return only 'peaks', 'valleys', or 'both'?
#' @param correctionRangeSeconds the range in which the real maximum/minimum value should be searched, in seconds.
#' around the derivative shift. Should be less then the periodicity of the signal.  0.5 is good for skin conductance.
#' @param minPeakDelta the minimum delta from valley to peak for detection. Skin conductance traditionally uses 0.05uS, or 0.03uS
#'
#' @return a list of: "bool" a logical vector of the same length of x with TRUE value corresponding to a match;
#' "samples" the index of matches in x;
#' "seconds" the position in seconds of matches (if x is a ts object);
#' and "type" a charachter vector defining for each "samples" value if its a 'p' peak or a 'v' valley.
#' @export
peakFinder = function(x, sgol_p = 2, sgol_n = 25, mode=c("peaks","valleys","both"), correctionRangeSeconds, minPeakDelta){
  sampRate = frequency(x)
  if(is.na(sgol_p)||is.na(sgol_n)){
    smooth_x = x
  } else {
  smooth_x = signal::sgolayfilt(x,  p =sgol_p, n = sgol_n, m=0)
  }
  fdderiv1  = diff(smooth_x)
  
  mode = match.arg(mode)
  pik = sign(embed(fdderiv1,2)) #embed appaia al segnale il segnale laggato di 1 samp
  s = pik[,2] - pik[,1] #that's where the magic happens
  s[is.na(s)] = 0
  
  #always both
  pikboo = c(FALSE,abs(s) ==  2, FALSE)
  piksam = which(pikboo) 
  s = s[s!=0]
  pv = s
  pv[s>0] = "p"
  pv[s<0] = "v"
  #correzione manuale: cerca il valore pi- alto nei dintorni del cambio di derivata
  
  
  
  any(diff(piksam)<= 0)
  correctionRangeSamp = correctionRangeSeconds*sampRate
  for(v in seq_along(piksam)){


    #individua il range con 0.5s prima e dopo la valle della derivata (esclusi gli estremi inzio e fine ts)
    prevv = if(v==1) 1 else piksam[v-1] +1
    nextv = if(v==length(piksam)) length(x) else piksam[v+1] -1
    search_interval = 
      seq(from = max(1,piksam[v]-correctionRangeSamp, prevv),
          to   = min(length(x),piksam[v]+correctionRangeSamp, nextv),
          by=1)
    
    if(FALSE) {
      plot(piksam[v]+((-50):(50)),x[piksam[v]+((-50):(50))], main=v)
      points(piksam, x[piksam],col=2,pch=18)
      text(piksam, eda[piksam]+0.2, paste(1:length(piksam), pv),cex=0.7)
      points(piksam[v],x[piksam[v]],col=3,pch=4,cex=3)
      points(range(search_interval),x[range(search_interval)],col=4,pch=4,cex=3)
      
    }
 
    #always both
    piksam[v] = search_interval[which.minmax(x[search_interval],pv[v])]
    piksam[v] = max(1, piksam[v])

    piksam[v] = min(piksam[v], length(x))
    if(FALSE){
      points(piksam[v],x[piksam[v]],col="gold",pch=17,cex=2)
    }
  }

  #trova picchi con sd pi- bassa di tot, sono solo rumore in un segnale essenzialmente piatto
  #idee: fisso a IQR(x)/20 ma magari si pu- trovare un valore pi- sensato
  #     -fisso a 0.05uS da onset a picco
  toDelete=c()
  for(v in 2:(length(piksam)-1)){
    if(pv[v]=="p"){
      search_interval = x[piksam[v-1]:piksam[v+1]]
      search_interval = search_interval[!is.na(search_interval)]
      # plot(c(10,16,search_interval),t="l",main=paste(v," - sd:",round(sd(search_interval,na.rm=T),2)))
      if(length(search_interval) == 0 || (max(search_interval)-search_interval[1])<minPeakDelta){#delta dall'onset al picco almeno 0.05uS o tutti NA
        # if(sd(search_interval)<IQR(x)/20){
        #se alcuni picchi vanno rimossi perch- sono essenzialmente flat
        #non ci possono essere 2 valley di seguito, quindi elimina quella col valore pi- alto
        
        fakeValley = c(v-1,v+1)[which.max(c(x[piksam[v-1]] ,x[piksam[v+1]]))]
        #elimina sia il picco 'v' che la fake valley
        toDelete = c(toDelete, v,fakeValley)
      }
    }
  }
  if(length(toDelete)>1){
    piksam = piksam[-toDelete]
    pv = pv[-toDelete]
    }
  
  
  #se ci sono tante valley, togli intanto tutte le v che hanno v sia a destra che sinistra
  toDelete=c()
  for(v in 2:(length(pv)-1)){
    # lel = pv=="v"
    # lele = embed(lel,2)
    # err = which(lele[,2] + lele[,1]>1)
    if(pv[v]=="v" && pv[v-1]=="v" && pv[v+1]=="v"){
      toDelete = c(toDelete, v)
    }
  }
  if(length(toDelete)>1){
    piksam = piksam[-toDelete]
    pv = pv[-toDelete]}
  
  #ora ci saranno al massimo dei casi p--v--v--p
  #eliminale solo se poco distanti
  toDelete = c()
  for(v in 1:(length(pv)-1)){
    if(pv[v]=="v" && pv[v+1] =="v"){
      # print((piksam[v+1]-piksam[v])/sampRate)
      # if ((piksam[v+1]-piksam[v]) < 5*sampRate  ){
        toDelete = c(toDelete, c(v+1,v)[which.max(c(x[piksam[v+1]],x[piksam[v]]))])
      # }
    }
  }
  if(length(toDelete)>1){
    piksam = piksam[-toDelete]
    pv = pv[-toDelete]}
  
  
  #tieni solo picchi e valli, se vuoi!
  if(mode=="peaks") {
    piksam = piksam[which(pv=="p")]
    pv = pv[which(pv=="p")]
  } else if(mode == "valleys"){
    piksam = piksam[which(pv=="v")]
    pv = pv[which(pv=="v")]
    
    
  }
  
  #ricrea pikboo & piks dai sample corretti
  pikboo = rep(F,length(pikboo))
  pikboo[piksam] = T
  piks = time(x)[piksam]
  
  list("bool" = pikboo,
       "samples" = piksam,
       "seconds" = piks,
       "type" = pv,
       "y" = x[piksam])
}
which.minmax = function(x, pv){if(pv=="p") which.max(x) else if(pv=="v") which.min(x)}

