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




#' Peak matching similarity of signals
#'
#' @param experiment 
#' @param signal the name of a DyadSignal contained in experiment 
#' @param lagSec 
#' @param weightMalus weightmalus è la percentuale di malus per il lag più estremo. Es, con weightMalus= 20 e r = 1 al massimo lag, la trasformazione diventa r' = 0.8
#' @param match_threshold a value between 0 and 1, specifying the similarity threshold to assign a peak-peak match.
#' @param minSizeSec the smallest allowed window size. very small windows may artificially inflate synchrony.
#' @param algorithm character. The version of AMICo algorithm. 
#' @param outputName
#' @param ... Further arguments passed to peakFinder. Specifically, 
#' sgol_p, sgol_n, correctionRangeSeconds and minPeakAmplitude must be set
#' according to the signal type.
#' @param maxPeakDuration 
#' @param interval_sec in AMICO2, how long should unmatched sequences be?
#'
#' @details Good values for skin conductance are correctionRangeSeconds = 0.5,
#' minPeakAmplitude = 0.05.
#' 
#' The difference in algorithms can be found in the changelog. Basically between
#' 1.0 and 1.1 there is a change in peakFinder to avoid detecting micro-peaks in generally flat areas.
#' 
#' The difference with 2.0 is much more substantial. Instead of calculating
#' similarities twice, in v2 the same approach is used to choose the best peak-to-peak
#' associations and to report the final values. Also exclusively valley-peak-valley
#' or valley-peak-halfrecovery sequences are used for the calculation.
#' Unmatched sequences are evenly split at interval_sec intervals.
#' 
#' XBEST guide:
#' row: the peak number of s1 wich was matched
#' col: the peak number of s2
#' similarity: some similarity computation
#' lag: the delta between the sample of p1 and p2
#' a1, p1, b1: sample position of onset, peak, and end of a feature
#' a2, p2, b2: the same for s2
#' a, b: the average start and end between s1 and s2
#' ta, tb: the time corresponding to a and b
#' 
#' @return
#' @export
#'
#' @examples
AMICo = function(experiment, signal="all", lagSec=6, #@OBSOLETE 2022 pmBest diventa AMICo
                  weightMalus = 50,
                  match_threshold = 0.1, minSizeSec=4,
                  maxPeakDuration = "hrt", interval_sec = 10,
                  algorithm=c("v1.1", "v1.0", "v2RC0.4"),
                  outputName = "AMICo",
                  sgol_p=2,sgol_n=25,correctionRangeSeconds,minPeakAmplitude){
  
  #debug
  # experiment = lr; signal="SC"; lagSec=4;
  # sgol_p = 2; sgol_n = 25;  weightMalus = 35;
  # match_threshold = 0.25; minSizeSec=1;
  # algorithm=c("AMICo_v1.0");
  # outputName = "PMBest";
  # correctionRangeSeconds = 0.5; minPeakAmplitude = 0.05;
  


  algorithm=match.arg(algorithm)
  # if(grepl("v2",algorithm)){
  #   warning("V2 is not yet fully implemented. some hidden paramenters are default")
  # }
  if(!is(experiment,"DyadExperiment")) stop("Only objects of class DyadExperiment can be processed by this function")

  nSessions = length(experiment)

  # progresbar
  pb <- progress::progress_bar$new(
    format = "Calculation::percent [:bar] :elapsed | ETA: :eta",
    total = nSessions,    # number of iterations
    width = 60, 
    show_after=0 #show immediately
    )
  
  progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar
  
  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  } 
  
  opts <- list(progress = progress)
  
  
  #parallelization
  warningsFile = "AMICoWarnings"
  if (file.exists(warningsFile)) {
    unlink(warningsFile)
  }
  cores=parallel::detectCores()-1
  cat(paste0("\r\nPerforming parallelized computation of AMICo ",algorithm," using ",cores," cores.\r\n")) #verified!
  cat(paste0("\r\nHigh Sync at positive lags implies that the ",s2Name(experiment[[1]]),
             " follows the ", s1Name(experiment[[1]]),"\r\n")) #verified!
  cl <- parallel::makeCluster(cores[1], outfile=warningsFile) #not to overload your computer
  doSNOW::registerDoSNOW(cl)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  pb$tick(0)
  #Questo è un loop parallelo sulle sessions.
  #il return del ciclo deve essere un oggetto DyadSession
  experiment2 <- foreach::foreach(
    iSession = 1:nSessions,
    .options.snow = opts, .errorhandling='pass'
    ) %dopar% { 
      
      xsignal = experiment[[iSession]][[signal]]
      
      if(grepl("v2",algorithm)){
        #I am using v2 which already implements peakmatch and ppsync in the same
        #function

        synchro = AMICo2(xsignal, lagSec=lagSec, 
                            weightMalus=weightMalus, match_threshold=match_threshold,
                        minSizeSec=minSizeSec,interval_sec=interval_sec,
                            outputName=outputName, maxPeakDuration=maxPeakDuration,
                        sgol_p=sgol_p,sgol_n=sgol_n,
                        correctionRangeSeconds=correctionRangeSeconds,
                        minPeakAmplitude=minPeakAmplitude )
      } else {
        synchro = AMICo1(xsignal, lagSec=lagSec, 
                            weightMalus=weightMalus, match_threshold=match_threshold,
                            outputName=outputName,
                            algorithm=algorithm, minSizeSec=minSizeSec,
                         sgol_p=sgol_p,sgol_n=sgol_n,
                         correctionRangeSeconds=correctionRangeSeconds,
                         minPeakAmplitude=minPeakAmplitude)
        
      }
      experiment[[iSession]][[signal]][[outputName]] = synchro
      return(experiment[[iSession]])
    }
  parallel::stopCluster(cl) 
  
  #restore session names and attributes
  for(iSession in 1:nSessions){
    if(!is.null(experiment2[[iSession]][["message"]])){
      write(paste0("QQQZ",experiment2[[iSession]][["message"]],"\n"),file=warningsFile,append=TRUE)
    }
    experiment2[[iSession]] = cloneAttr(experiment[[iSession]],experiment2[[iSession]])
  }


  cat("\r\nDone ;)\r\n")
  names(experiment2) = names(experiment)
  experiment2 = cloneAttr(experiment, experiment2)
  
  
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
AMICo1 = function(signal,lagSec,weightMalus,match_threshold, outputName, algorithm,
                  minSizeSec,...) {
  args <- c(as.list(environment()), list(...),
            list(
              sessionId=sessionId(signal),
              dyadId=dyadId(signal),
              groupId=groupId(signal))
  )
  args["signal"] = NULL
  l = list(...)
  sgol_p=l[["sgol_p"]]
  sgol_n=l[["sgol_n"]]
  correctionRangeSeconds=l[["correctionRangeSeconds"]]
  minPeakAmplitude = l[["minPeakAmplitude"]]
  
  
  xbest = peakMatch(signal, lagSec=lagSec, 
                      weightMalus=weightMalus, match_threshold=match_threshold,
                      outputName=outputName,algorithm=algorithm, ...)
  #the second calculates the actual sync values
  res = ppSync_dev(signal,xbest,minSizeSec,outputName=outputName)
  
  return(newAMICo(res$sync, res$lag, res$xBest, args))
  
}


## peak picking best lag
peakMatch = function(signal,lagSec=4, sgol_p = 2, sgol_n = 25, weightMalus = 30,
                     match_threshold=0.0, outputName, correctionRangeSeconds,
                     minPeakAmplitude, ...) {

  if(!is(signal,"DyadSignal")) stop("Only objects of class DyadSignal can be processed by this function")
  
  #signal = lr$all_vid1_01$contempt
  # signal = lr10$all_CC_1$SC
  # signal = lr10f$all_02_1$SC
  #pacs best settings
  ## lr = pmBest(d,signals = "all",lagSec=4,match_threshold=0.25,
  #             minSizeSec=5,weightMalus = 35 ,algorithm = "AMICo",outputName = "PMdev")
  
  d = signal$s1
  d2  = signal$s2
  SR = frequency(signal)
  ransamp = lagSec * SR
  
  
  ### peaks-valleys detection ########
  l = list(...)
  if(any(names(l) == "algorithm") && l[["algorithm"]]=="v1.0"){
    print("v1.0 legacy")
    s1b = legacyPeakFinder(d,  sgol_p, sgol_n, mode = "b", correctionRangeSeconds)
    s2b = legacyPeakFinder(d2, sgol_p, sgol_n, mode = "b", correctionRangeSeconds)
  } else {
    s1b = peakFinder(d,  sgol_p, sgol_n, mode = "b", correctionRangeSeconds, minPeakAmplitude)
    s2b = peakFinder(d2, sgol_p, sgol_n, mode = "b", correctionRangeSeconds, minPeakAmplitude)
  }
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
  # signal[[outputName]] = newAMICo(NULL,NULL,xbest, NULL)
  # signal
  xbest
}


#AMICo: Adaptive Matching Interpolated Correlation
ppSync_dev = function(signal, xbest, minSizeSec, outputName) {
  ##questa è una versione di debug e sviluppo di ppSync di rIP
  # richiede dei dati generati con rIP in un oggetto chiamato lr
  # ottimi quelli di "PACS_rIP engine.R"
  #
  #la versione attuale (l'ultima su github ) è imperfetta:
  # 1. usare la finestra più lunga causa degli errori in caso di cambio di lag:
  #   ^  ^
  # /  \/ \
  # __^.^__    <-qui anche se i due picchi sono uguali la cor viene bassa
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
  # if(is.null(signal[[outputName]])||is.null(signal[[outputName]]$xBest))stop("Please run peakMatch before.")
  if(nrow(xbest)<2) stop("peakMatch found only one peak in the signal")
  if(nrow(xbest)<10) warning("peakMatch found less than 10 peaks")
  d = signal$s1
  d2  = signal$s2
  lagSec = attr(signal[[outputName]], "lagSec")
  SR = frequency(signal)
  ransamp = lagSec * SR
  # xbest = signal[[outputName]]$xBest
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
    
    if(length(rs1) > minSizeSec*SR){
      # #normalizing element
      # if(scaledCorrelation){
      #   # #his is shit
      #   # rs1=c(nMin,nMax,rs1)
      #   # rs2=c(nMin,nMax,rs2)
      #   rs1 = rs1-mean(rs1)
      #   rs1 = rs1/normPam
      #   rs2 = rs2-mean(rs2)
      #   rs2 = rs2/normPam
      #   
      #   iCor = sqrt(sum((rs2-rs1)^2))
      #   iCor = 1- iCor/maxDist
      # }else 
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
  
  
  # finallyt instantiate new sync class object
  sync = rats(syncvec, start = start(d),
              frequency=frequency(signal), timeUnit="second",
              valueUnit="Pearson correlation")
  lags = rats(lagvec,  start = start(d),
              frequency=frequency(signal), timeUnit="second",
              valueUnit="seconds")
  # applica gli artefatti
  for(i in seq_len(nrow(signal$artefacts)) ){ 
    window(sync, signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
    window(lags, signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
  }
  
  return(
    list (
     xBest = xbest,
     sync = sync,
     lags = lags
     )
  )
  
}

#' Finds peaks and throughs features in a time series
#'
#' @param x a rats, ts, or a numeric vector.
#' @param sgol_p smoothing filter order. If NA, the filtering is disabled.
#' @param sgol_n smoothing filter length (must be odd).
#' @param mode should the function return only 'peaks', 'valleys', or 'both'?
#' @param correctionRangeSeconds the range in which the real maximum/minimum
#' value should be searched, in seconds.
#' around the derivative shift. Should be less then the periodicity of the
#' signal.  0.5 is good for skin conductance.
#' @param minPeakAmplitude the minimum delta from valley to peak for detection.
#' Skin conductance traditionally uses 0.05uS, or 0.03uS
#' @details The function has a relatively simple yet robust implementation:
#' after a Sgolay smoothing, the first derivative's sign changes are detected.
#' Then, the actual minima or maxima is looked for in a given range.
#' NOTE: the function ensures that no consecutive throughs or peaks are returned.
#' 
#'
#' @return a list of:
#' \describe{
#'   \item{bool}{a logical vector of the same length of x with TRUE
#'    value corresponding to a match
#'    }
#'   \item{samples}{the position of features relatively to x's index}
#'   \item{time}{temporal coordinates of the features (if x has a frequency attribute)}
#'   \item{type}{a character vector defining for each "samples" value if its a 'p' peak or a 'v' valley (trough)}
#'   \item{y}{The value of x correspoding to the detected features}
#'   \item{amp}{the positive (v->p) and negative (p->v) amplitudes of the detected features}
#'   \item{i}{The absolute position of the features along the signal. The values
#'      are independent from the 'mode' argument, allowing interaction between
#'      different calls of the function.}
#'   
#'   
#' }
#' @export
peakFinder = function(x, sgol_p = 2, sgol_n = 25, mode=c("peaks","valleys","both"),
                      correctionRangeSeconds, minPeakAmplitude){
  ####debug
  # stop("debug")
  # x = d1
  # sgol_p = 2
  # sgol_n = 25
  # mode="both"
  # correctionRangeSeconds = 0.5
  # minPeakAmplitude = 0.05
  ##
  
  if(missing(sgol_p)) stop("in peakfinder missing sgol_p")
  if(missing(sgol_n)) stop("in peakfinder missing sgol_n")
  if(missing(correctionRangeSeconds)) stop("in peakfinder missing correctionRangeSeconds")
  if(missing(minPeakAmplitude)) stop("in peakfinder missing minPeakAmplitude")
  
  SR = frequency(x)
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
  s = s[abs(s)== 2]
  pv = s
  pv[s>0] = "p"
  pv[s<0] = "v"
  #correzione manuale: cerca il valore pi- alto nei dintorni del cambio di derivata
  
  
  
  # any(diff(piksam)<= 0)
  correctionRangeSamp = correctionRangeSeconds*SR
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
      if(length(search_interval) == 0 || (max(search_interval)-search_interval[1])<minPeakAmplitude){#delta dall'onset al picco almeno 0.05uS o tutti NA
        # if(sd(search_interval)<IQR(x)/20){
        #se alcuni picchi vanno rimossi perch- sono essenzialmente flat
        #non ci possono essere 2 valley di seguito, quindi elimina quella col valore pi- alto
        
        fakeValley = c(v-1,v+1)[which.max(c(x[piksam[v-1]] ,x[piksam[v+1]]))]
        #elimina sia il picco 'v' che la fake valley
        toDelete = c(toDelete, v,fakeValley)
      }
    }
  }
  if(length(toDelete)>0){
    piksam = piksam[-toDelete]
    pv = pv[-toDelete]
  }
  
  #risolvi il problema delle valli consevutive v-v-v...
  #this is ultrafast but don't discriminate vs and ps  
  # toDelete = which(pv[-length(pv)] == pv[-1])
  toDelete=c()
  for(v in 2:(length(pv))){
    #se la feature attuale è uguale alla precedente
    if(pv[v] == pv[v-1]){
      if(pv[v] == "v" ){
        #se la feature è valley, tieni solo l'ultima
        #ovvero elimina quella precedente
        toDelete = c(toDelete, v-1)
      } else if(pv[v] == "p"){
        #se la feature è un picco, tieni il più alto
        #ovvero elimina il più basso.
        #se il più basso è v, lowest = 0, se no lowest = 1
        lowest = which.min(c(x[v], x[v-1]))
        toDelete = c(toDelete, v-lowest)
      }
    }
  }
  if(length(toDelete)>0){
    piksam = piksam[-toDelete]
    pv = pv[-toDelete]
  }
  
  is = 1:length(piksam)
  
  #tieni solo picchi e valli, se vuoi!
  if(mode=="peaks") {
    piksam = piksam[which(pv=="p")]
    is = is[which(pv=="p")]
    pv = pv[which(pv=="p")]
    
  } else if(mode == "valleys"){
    piksam = piksam[which(pv=="v")]
    is = is[which(pv=="v")]
    pv = pv[which(pv=="v")]
  }
  
  #ricrea pikboo & piks dai sample corretti
  pikboo = rep(F,length(pikboo))
  pikboo[piksam] = T
  piks = time(x)[piksam]
  
  
  list("bool" = pikboo,
       "samples" = piksam,
       "time" = piks,
       "class" = pv,
       "y" = x[piksam],
       "amp" = c(NA,diff(x[piksam])),
       "index" = is
  )
  
}
which.minmax = function(x, pv){if(pv=="p") which.max(x) else if(pv=="v") which.min(x)}




