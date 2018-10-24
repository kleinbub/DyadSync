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
#' @return
#' @export
#'
#' @examples
pmBest = function(experiment, signals="all", lagSec=7,
                  sgol_p = 2, sgol_n = 25,  weightMalus = 30,
                  match_threshold = 0.0,minSizeSec=5, algorithm=c("classic","dev","sccf"), outputName = "PMBest"){
  algorithm=match.arg(algorithm)
  if(!is(experiment,"DyadExperiment")) stop("Only objects of class DyadExperiment can be processed by this function")
  cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]),
             " follows the ", s1Name(experiment[[1]]),"\r\n")) #verified!
  nSessions = length(experiment)
  experiment2 = Map(function(session,iSession){
    if(signals=="all") signals = names(session)
    cat("\r\n",paste(dyadId(session),session(session)))
    session[signals] = Map(function(signal){
      cat(" |",name(signal))
      #the first function calculates best lag
      signal = peakMatch(signal, lagSec=lagSec, sgol_p=sgol_p, sgol_n=sgol_n,
                         weightMalus=weightMalus, match_threshold=match_threshold, outputName=outputName)
      #the second calculates the actual sync values
      if(algorithm == "dev")
        signal = ppSync_dev(signal,minSizeSec,outputName=outputName)
      else if(algorithm == "classic")
        signal = ppSync(signal,minSizeSec,outputName=outputName)
      else if(algorithm == "sccf")
        signal = ppSync_sccf(signal,outputName=outputName)
      
      return(signal)
    }, session[signals])
    return(session)
  },experiment,seq_along(experiment))
  cat("\r\nDone ;)\r\n")
  classAttr(experiment2)=classAttr(experiment)
  return(experiment2)
}


## peak picking best lag
peakMatch = function(signal,lagSec=8, sgol_p = 2, sgol_n = 25, weightMalus = 30,match_threshold=0.0, outputName) {
  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  
  # signal = loncc$CC_1$SC
  
  d = signal$s1
  d2  = signal$s2
  sampRate = frequency(signal)
  ransamp = lagSec * sampRate
  
  sd1 = c(0,rangeRescale(diff(d),-0.5,0.5,-max(abs(diff(d))),max(abs(diff(d))) ))
  sd2 = c(0,rangeRescale(diff(d2),-0.5,0.5,-max(abs(diff(d2))),max(abs(diff(d2))) ))
  
  ### peaks-valleys detection ########
  s1p = peakFinder(d,  sgol_p, sgol_n, mode = "p", 0.5)
  s1v = peakFinder(d,  sgol_p, sgol_n, mode = "v", 0.5)
  s2p = peakFinder(d2, sgol_p, sgol_n, mode = "p", 0.5)
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
      for(f in 1:nMatch){
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
  xbest$lag =   allpik2[xbest$col]-allpik1[xbest$row]
  xbest$s1 = allpik1[xbest$row]
  xbest$s2 = allpik2[xbest$col]
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
    ii = i
    ab = list()
    ab[[1]] = xbest$s1[i]:xbest$s1[i+1]
    ab[[2]] = xbest$s2[i]:xbest$s2[i+1]
    
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
  signal[[outputName]]$xBest = xbest
  signal[[outputName]]$sync = DyadStream(syncvec, "PMBest_Sync", col=rgb(191,50,59,max=255))
  signal[[outputName]]$lag = DyadStream(lagvec, "PMBest_Lag", col=rgb(253,177,2,max=255))
  return(signal)
  
}

ppSync_dev = function(signal,minSizeSec, outputName) {
  ##questa è una versione di debug e sviluppo di ppSync di rIP
  # richiede dei dati generati con rIP in un oggetto chiamato lr
  # ottimi quelli di "PACS_rIP engine.R"
  #
  #la versione attuale (l'ultima su github ) è imperfetta:
  # 1. usare la finestra più lunga causa degli errori in caso di cambio di lag:
  #  ^  ^
  # / \/ \
  #   ^^    <-qui anche se i due picchi sono uguali la cor viene bassa
  # meglio se ALMENO una finestra sia sufficientemente lunga e quella breve venga interpolata
  #
  # 2. aggregando finestre diverse per avere una lunghezza minima si sminchia l'appaiamento fra picchi
  # cosa che si vede già con ii = c(2,3) {II e III righe di xbest} che vengono aggregate ma il lag risultante è errato.
  
  ## idee: - pesare per la differenza di durata?
  ##       - pesare per la differenza di ampiezza (normalizzata sulla SD individuale?)
  
  cat(" - dev version of ppSync")
  #please refer to ppSync dev.R in research folder of DyadClass
  
  # signal = lr$CC_1$SC
  # minSizeSec=5 <-- can be smaller as the check is on the shorter segment?
  if(is.null(signal[[outputName]])||is.null(signal[[outputName]]$xBest))stop("Please run peakMatch before.")
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
      ab = (xbest$s1[i-1]+1):xbest$s1[i]
      lagvec[ab] = round(seq(from=lags[i-1], to=lags[i], length.out = length(ab)))
    }
  }
  lagvec[is.na(lagvec)] = lags[length(lags)] #fix the last values
  
  
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
    
    if(length(rs1) >= minSizeSec*sampRate){
      iCor = cor(rs1,rs2, use = "c")
      # iCor = cor(diff(rs1),diff(rs2), use = "c")
      
      # ### diagnostic plot ####
      # par(mfrow=c(2,2))
      # 
      # # what the correlation sees:
      # rs1p = rangeRescale(rs1,0,1)
      # rs2p = rangeRescale(rs2,0,1)
      # plot(rs1p,type="l",main=iCor,ylim=c(0,1), col=attr(signal$s1,"col"))
      # lines(rs2p, col=attr(signal$s2,"col"))
      # ##what I see:
      # s1range = xbest$s1[iCum[1]]:xbest$s1[iCum[length(iCum)]+1]
      # s2range = xbest$s2[iCum[1]]:xbest$s2[iCum[length(iCum)]+1]
      # s1y= d[s1range]; s2y = d2[s2range]
      # plot(s1range,s1y,
      #      xlim=c(min(s1range,s2range),max(s1range,s2range)), ylim=c(min(s1y,s2y),max(s1y,s2y)),
      #      type="l",main=iCor, col=attr(signal$s1,"col"))
      # lines(s2range,s2y, col=attr(signal$s2,"col"))
      # 
      # ##dividing by whole signal sd
      # plot(rs1x,type="l",main=cor(rs1x,rs2x),ylim=c(0,1), col=attr(signal$s1,"col"))
      # lines(rs2x, col=attr(signal$s2,"col"))
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
  signal[[outputName]]$sync = DyadStream(lcc, "PMBest_Sync", col=rgb(191,50,59,max=255), start=xStart,frequency=1/incSec )
  signal[[outputName]]$lag = DyadStream(lagvec, "PMBest_Lag", col=rgb(253,177,2,max=255), start=start(signal),frequency=sampRate )
  return(signal)
  
}

##peakfinder è una funzione ausiliaria, è sviluppata in "peak-centered-flex-CCF_v1.93"
#' Title
#'
#' @param x a ts object or a numeric vector.
#' @param sgol_p smoothing filter order.
#' @param sgol_n smoothing filter length (must be odd).
#' @param mode should the function return only 'peaks', 'valleys', or 'both'?
#' @param correctionRangeSeconds the range in which the real maximum/minimum value should be searched, in seconds.
#' around the derivative shift. Should be less then the periodicity of the signal.  0.5 is good for skin conductance.
#' @param valid a logical vector of the same length of x. No peaks/valleys are found where valid is FALSE.
#'
#' @return a list of: "bool" a logical vector of the same length of x with TRUE value corresponding to a match;
#' "samples" the index of matches in x;
#' "seconds" the position in seconds of matches (if x is a ts object);
#' and "type" a charachter vector defining for each "samples" value if its a 'p' peak or a 'v' valley.
#' @export
peakFinder = function(x, sgol_p = 2, sgol_n = 25, mode=c("peaks","valleys","both"), correctionRangeSeconds = 0.5, valid){
  if(missing(valid)) valid = rep(T,length(x))
  sampRate = frequency(x)
  smooth_x = signal::sgolayfilt(x,  p =sgol_p, n = sgol_n, m=0)
  fdderiv1  = diff(smooth_x)
  mode = match.arg(mode)
  pik = sign(embed(fdderiv1,2)) #embed appaia al segnale il segnale laggato di 1 samp
  s = pik[,2] - pik[,1] #that's where the magic happens
  if(mode=="peaks")
    pikboo = c(s ==  2, FALSE) #embed perde 1 sample, quindi aggiungi un FALSE alla fine
  else if(mode=="valleys")
    pikboo = c(s == -2, FALSE) #embed perde 1 sample, quindi aggiungi un FALSE alla fine
  else pikboo = c(abs(s) ==  2, FALSE)
  
  pikboo[valid==F] = FALSE #cancella i picchi nelle aree con artefatti
  
  piksam = which(pikboo) #in quali sample c'è un'inversione di segno della derivata?
  #correzione manuale: cerca il valore più alto nei dintorni del cambio di derivata
  if(mode=="peaks")
    pv = rep("p",length(piksam))
  else if(mode=="valleys")
    pv = rep("v",length(piksam))
  else {
    s = s[s!=0]
    pv = s
    pv[s>0] = "p"
    pv[s<0] = "v"
  }
  
  for(v in seq_along(piksam)){
    #individua il range con 0.5s prima e dopo la valle della derivata (esclusi gli estremi inzio e fine ts)
    search_interval = x[
      max(1,piksam[v]-correctionRangeSeconds*sampRate):min(length(x),piksam[v]+correctionRangeSeconds*sampRate)
      ]
    #trova il più piccolo e aggiorna pikmsam
    if(mode=="peaks")
      piksam[v] = max(1, piksam[v]+  (which.max(search_interval) - round(length(search_interval)/2)))
    else if (mode=="valleys")
      piksam[v] = max(1, piksam[v]+  (which.min(search_interval) - round(length(search_interval)/2)))
    else
      piksam[v] = max(1, piksam[v]+  (which.minmax(search_interval,pv[v]) - round(length(search_interval)/2)))
    
    piksam[v] = min(piksam[v], length(x))
  }
  #ricrea pikmboo dai sample corretti
  pikboo = rep(F,length(pikboo))
  pikboo[piksam] = T
  piks = time(x)[piksam]
  list("bool" = pikboo,
       "samples" = piksam,
       "seconds" = piks,
       "type" = pv)
}
which.minmax = function(x, pv){if(pv=="p") which.max(x) else if(pv=="v") which.min(x)}

