# Dude I have checked this function like, 100 times, it is very straightforward.
# there is no issue in it. Go search somewhere else.
# @OBSOLETE forceMaxPads parameter

#' @param x a DyadExperiment Object.
#'
#' @param plotPrefix character. path and prefix of output plots names. If NA, no plot will be drawn.
#' @param signal character. name of a DyadSignal contained in x.
#' @param sync character. Name of a synchronization object. E.g. PMBest, or CCFBest
#' @param series character. Name of a rats object contained within sync
#' @param category character. Name of a DyadCategory object contained in x
#' @param categoryIndex character. Name of a factor column in 'category'.
#' @param absolute logical. should the series be transformed to absolute values?
#' @param minEpochSec number of seconds of minimum epoch. Shorter epochs will be deleted.
#' @param minOccurrences min occurrences of given category to be kept. Less frequent categories will be deleted.
#' @param keepOnly vector of characters denoting if only given categories must be analyzed. Eg: c("cat1","cat2").
#' @param nIter integer. Number of random extractions.
#' @param rotate if true windows are extracted from a randomly shifted timeseries. This is better TRUE.
#' @param extract_from Either "all", or "remaining". Should new extractions be done randomly from the whole series or only on non-categorized parts.
#' @param epochSummaryFUN String. the function used to summarize EACH epoch's time-series to a single value. Suggested: median.
#' @param centralityFUN  String. A centrality function used to synthetize the distribution of ALL epochs. Suggested: median.
#' @param dispersionFUN  String. A dispersion function used to synthetize the distribution of ALL epochs. Suggested: MAD
#' @param forceMaxPads logical. If true the random epochs are evenly spaced across the whole data. This is for retrocompatibility only. Keep disabled.
#' @param returnData logical. Should the data be returned?
#' @param credibilityMass 
#' @param nLines 
#' @param xlim either "auto" or a numeric vector of length 2 to set the plot scale and density boundaries.
#' @param ... 
#' 
#' @details Credible intervals are calculated with \link[hdr]{hdrcde} function
#' @return if returnData is TRUE, a list with all observed and random data and statistics.
#' If plotPrefix is not NA, a plot will be drawn too.
#' @export
#'
#' @examples
categoryPerm = function(x,
                        plotPrefix,
                        signal,
                        sync,
                        series,
                        category,
                        categoryIndex,
                        keepOnly = c(), #
                        minEpochSec, #
                        minOccurrences, #
                        
                        absolute = FALSE,
                        returnData = FALSE,
                        
                        #randomization settings
                        nIter = 1000,
                        rotate = TRUE, 
                        extract_from = c("all", "remaining")[2], 
                        epochSummaryFUN = c("mean","median")[2],
                        centralityFUN = c("mean","median")[2],
                        dispersionFUN = c("sd","mad")[2],
                        forceMaxPads = FALSE,
                        
                        #analyses settings
                        credibilityMass = 0.89,
                        
                        #graphic settings
                        nLines=200,
                        xlim = "auto",
                        ...) {
  
  args = c(as.list(environment()), list(...))
  args$x = NULL
  ####DEBUG
  # rm(list=ls())
  # load("A:/OneDrive - Università degli Studi di Padova/__Ricerche/2022_IM_paper/IM_engine_2023_amico2_v2.RData")
  # plotPrefix = paste0("permplottest_")
  # signal="SC"
  # sync="amico2";sync="PMdev"
  # series="sync"
  # category = "IM"
  # categoryIndex = "tipo"
  # # keepOnly = c("SI") #
  # minEpochSec = 3 #
  # minOccurrences = 5 #
  # returnData = TRUE
  # absolute = FALSE
  # credibilityMass = 0.89
  # stop("debug")
  # #randomization settings
  # nIter = 100
  # rotate = TRUE
  # extract_from = c( "remaining")
  # epochSummaryFUN = c("median")
  # centralityFUN = "median"
  # dispersionFUN = "mad"
  # #graphic settings
  # xlim = c(-1,1);xlim = "auto"
  # d3=d2
  # keepOnly = c("Reconceptualization")
  ###
  
  d3 = x
  if(absolute){
    for(i in 1:length(d3)){
      d3[[i]][[signal]][[sync]][[series]] = abs(d3[[i]][[signal]][[sync]][[series]])
    }
  }
  samplesPerSec = sapply(d3, function(x){frequency(x[[signal]][[sync]][[series]])})
  samplesPerSec = unique(samplesPerSec)
  if(length(samplesPerSec)>1) stop("multiple frequencies detected for selected series in different DyadSessions:\r\n",samplesPerSec)
  cat("\r\nSEGMENTING BY EPOCH:\r\n")
  d4 = epochSeries(d3, signal=signal, sync=sync,series = series,
                   category=category,categoryIndex=categoryIndex, mergeEpochs = F, artefact.rm=F)
  
  resName = paste0(c(toupper(category), "_", totitle(c(sync,series))),collapse = "")
  ex2 = extractEpochs (d4, signal=signal, sync=sync, series=series,
                       category=category, categoryIndex=categoryIndex)
  lN   = unlist(lapply(ex2,length))
  
  ## keep only selected categories
  if(length(keepOnly)>0){
    keepOnly = c(keepOnly, "remaining")
    ex2 = ex2[keepOnly]
  }
  ex2 = ex2[sapply(ex2,length)>0]
  

  
  shift <- function(x, n) {
    n = n%%length(x)
    if(n==0) x
    else c(x[(n+1):length(x)], x[1:n])
  }
  
  fastcohen = function(a,b){
    #this requires that a and b don't have missing values
    (mean(a)-mean(b))/sqrt(((length(a)-1)*sd(a)^2 +  (length(b)-1)*sd(b)^2 )/(length(a)+ length(b)-2 ) )
  }
  

  p_summaryFUN     = eval(parse(text=epochSummaryFUN))
  p_centralityFUN  = eval(parse(text=centralityFUN))
  p_dispersionFUN  = eval(parse(text=dispersionFUN))
  
  
  # na.omit.safe = function
  #remove epochs shorter than minEpochSec seconds
  
  toRemove = c()
  infoRm = numeric(length(ex2))
  infoNA = numeric(length(ex2))
  for(i in 1:length(ex2)){
    #remove NAs from the sequences, so that they don't influence results
    for(j in 1:length(ex2[[i]])){
      if(all(is.na(ex2[[i]][[j]]))) infoNA[i] = infoNA[i] + 1
      ex2[[i]][[j]] = ex2[[i]][[j]][!is.na(ex2[[i]][[j]])]
    }
    #after removal of NAs, check if remaining size is above threshold
    durations = as.numeric(sapply(ex2[[i]],length))
    toKeep = durations > (minEpochSec*samplesPerSec) #e.g. 5 seconds x 10 samplepersecond
    infoRm[i] = sum(!toKeep) - infoNA[i]
    ex2[[i]] = ex2[[i]][toKeep]
    
    #if any category has less than 10 instances, add to remove list
    if(length(ex2[[i]])<minOccurrences) toRemove = c(toRemove, i)
    
  }
  
  #remaining is a special category and should never be removed
  remi = which(names(ex2)=="remaining")
  if(remi %in% toRemove) toRemove = toRemove[toRemove!= remi]
  
  
  #create a nice report
  lLab = names(ex2)
  lN = lN[lLab]
  lN2   = unlist(lapply(ex2,length))
  lkept = !logical(length(ex2))
  lkept[toRemove] = FALSE
  report = data.frame("Index levels" = lLab, "original n"=lN,
             "all NAs" = infoNA, "too short"=infoRm, "final n" = lN2, "enough.obs"=lkept,
             row.names = NULL)
  report = report[!report[,1]=="remaining",]
  print(report,row.names = FALSE)
  #remove from removelist
  if(length(toRemove)>0) ex2 = ex2[-toRemove]
  
  toPlot = names(ex2)[names(ex2)!="remaining"]
  # toPlot = names(ex2)
  
  # #experimental shite
  # for(i in 1:length(ex2)){
  #   #remove NAs from the sequences, so that they don't influence results
  #   for(j in 1:length(ex2[[i]])){
  #     plot(ex2[[i]][[j]],t="l", main=names(ex2)[i],ylim=c(-1,1))
  # 
  #   }
  # }
  ################

  
  exRan = vector("list",length(toPlot))
  names(exRan) = toPlot
  cat("\r\nANALYZING:\r\n")
  for(tipo2 in toPlot){
    # tipo2 = toPlot[3] ## <-----------------------------------------------####
    # print(tipo2)
    #remove NAs should be superfloous
    ex2[[tipo2]] = ex2[[tipo2]][sapply(ex2[[tipo2]], function(x){!all(is.na(x))})] #da FALSE solo se tutta la finestra è NA
    #n seq di tipo CODIFICA
    n = length(ex2[[tipo2]])
    
    for(i in 1:n){
      #elimina valori mancanti
      ex2[[tipo2]][[i]] = ex2[[tipo2]][[i]][!is.na(ex2[[tipo2]][[i]])]
    }
    #elimina occorrenze se togliendo gli NA siamo più corti di minEpochSec
    ex2[[tipo2]] = ex2[[tipo2]][unlist(lapply(ex2[[tipo2]],length)) > minEpochSec*samplesPerSec]
    
    

    
    #1 crea sequenze di finestre consecutive dalle duration
    durations = lapply(ex2, function(x)sapply(x,length))
    durReal = durations[[tipo2]]
    
    #calcola mediana su ciascuna finestra dei dati reali
    obsEpochSummary = as.numeric(unlist(lapply(ex2[[tipo2]], p_summaryFUN)))
    obsEpochSummary = obsEpochSummary[!is.na(obsEpochSummary)]
    
    #calcola statistica centrale sui dati reali
    obsCentrality = p_centralityFUN(obsEpochSummary)
    obsDispersion = p_dispersionFUN(obsEpochSummary)
    
    if(extract_from =="all") {
      #estrai tutti i segnali sync delle sedute e togli i NA
      full_from = lapply(d4, function(x) as.numeric(x[[signal]][[sync]][[series]]))
      full_from = do.call("c",full_from)
      full_from = full_from[!is.na(full_from)]
    } else if(extract_from =="remaining")  {
      #usa i dati 'remaining'
      if(!"remaining"%in%names(ex2)) stop("no data in 'remaining' category, try to use extract_from = 'all'")
      full_from = do.call("c",ex2$remaining)
      full_from = full_from[!is.na(full_from)]
      dur_re = durations$remaining
      if(length(full_from)/samplesPerSec/60/60 <1) warning ("random data is being extracted from less than 1 hour of data. Consider using extract_from = 'all'")
      
    }
    max_from = length(full_from)
    
    
    ####### SETUP PARALLELIZATION
    # progresbar
    pb <- progress::progress_bar$new(
      format = "Calculation::percent [:bar] :elapsed | ETA: :eta",
      total = nIter,    # number of iterations
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
    warningsFile = "MyWarnings"
    if (file.exists(warningsFile)) {
      unlink(warningsFile)
    }
    cores=parallel::detectCores()-1
    cat(paste0("\r\nPerforming parallelized permutations of '",tipo2,"' categories using ",cores," cores.\r\n")) #verified!
    cl <- parallel::makeCluster(cores[1], outfile=warningsFile) #not to overload your computer
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    pb$tick(0)
    resList <- foreach::foreach(
      k = 1:nIter,
      .options.snow = opts, .errorhandling='pass'
    ) %dopar% {
      
      ################# HERE THE PARALLEL JOB
      
      #randomizza l'ordine delle durate
      durRand = sample(durReal)
      xstart = xend = c(1,numeric(length(durRand)-1))
      for(i in 2:n) {
        xstart[i] = xstart[i-1]+ durRand[[i-1]]
        xend[i-1] = xstart[i-1] +durRand[[i-1]] -1
      }
      xend[n] = xstart[n] + durRand[n]-1 #fix the last line
      
      #genera paddings casuali
      maxgap = max_from - sum(durRand)
      
      seed = runif(n,0,1)
      seed = seed/sum(seed) #all the seeds sum to 1
      pads = floor(seed*maxgap) #now all pads summed are circa maxgap
      if(!forceMaxPads){
        # @OBSOLETE this should just be the ONLY behaviour. remove the if.
        scaler = runif(n, 0,1) #each pad is scaled randomly from 0 to 100%
        pads = floor(pads*scaler)
      }

      
      #sposta xstart e xend
      rstart = xstart; rend =xend
      for(i in 1:n) {
        rstart[i:n] = rstart[i:n] + pads[i]
        rend[i:n] = rend[i:n] + pads[i]
      }

      durDelta = abs(sum(rend-rstart) - sum(durReal))
      if(durDelta > sum(durReal)*0.05) warning("random duration was very small")
      
      #ruota full_from a caso
      if(rotate) {
        rfrom_seed = round(runif(1,1,max_from))
        rfrom = shift(full_from, rfrom_seed)
      } else rfrom = full_from
      
      
      #crea un vettore di n finestre
      rfsplit = vector("list",n)
      for(i in 1:n){
        #seleziona da rfrom (ruotato) i valori delle finestre random
        rfsplit[[i]] = rfrom[(rstart[i]):(rend[i])]
        #elimina valori mancanti
        rfsplit[[i]] = rfsplit[[i]][!is.na(rfsplit[[i]])]
      }
      #elimina occorrenze se togliendo gli NA siamo più corti di minEpochSec
      rfsplit = rfsplit[unlist(lapply(rfsplit,length)) > minEpochSec*samplesPerSec]
      
      #estrai il parametro dalle finestre random
      #ora dovremmo avere la garanzia di avere dati senza NA, rendendo inutili gli na.rm=T
      randEpochSummary = as.numeric(unlist(lapply(rfsplit, p_summaryFUN)))
      randEpochSummary = randEpochSummary[!is.na(randEpochSummary)]
      
      randEpochSummary = as.numeric(unlist(lapply(rfsplit, p_summaryFUN)))
      randEpochSummary = randEpochSummary[!is.na(randEpochSummary)]
      
      if(length(xlim)==1 && xlim=="auto"){
        randDen = density(randEpochSummary)
      } else {
        randDen = density(randEpochSummary, from = xlim[1], to = xlim[2])
      }
      
      return(list(
        cohe = fastcohen(obsEpochSummary,randEpochSummary),
        cliff = effsize::cliff.delta(obsEpochSummary,randEpochSummary)$estimate,
        par1 = p_centralityFUN(randEpochSummary),
        par2 = p_dispersionFUN(randEpochSummary),
        randDen = randDen,
        randEpochSummary = randEpochSummary
        
      ))
      
    }
    parallel::stopCluster(cl) 
    
    #Check for issues
    for(i in 1:nIter){
      if(!is.null(resList[[i]][["message"]])){
        write(paste0("QQQZ",resList[[i]][["message"]],"\n"),file=warningsFile,append=TRUE)
      }
    }
    
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
    
    
    cat("\rExtracting values ...")
    #extract output elements
    cohenlist = unlist(lapply(resList, \(x){x$cohe}))
    cohenlist = cohenlist[!is.na(cohenlist)]
    
    clifflist = unlist(lapply(resList, \(x){x$cliff}))
    clifflist = clifflist[!is.na(clifflist)]
    
    ranCentraList  = unlist(lapply(resList, \(x){x$par1}))
    ranDisperList  = unlist(lapply(resList, \(x){x$par2}))
    
    randDenX = lapply(resList, \(x){x$randDen$x})
    randDenX = do.call("rbind",randDenX)
    randDenX = apply(randDenX,2,mean)
    
    minRandX = quantile(unlist(lapply(resList, \(x){min(x$randDen$x)})),0.05)
    maxRandX = max(unlist(lapply(resList, \(x){max(x$randDen$x)})),0.95)
    
    
    randDenY = lapply(resList, \(x){x$randDen$y})
    randDenY = do.call("rbind",randDenY)
    randDenY = apply(randDenY,2,mean)
    maxRandY = max(unlist(lapply(resList, \(x){max(x$randDen$y)})))
    
    pval =(length(ranCentraList[ranCentraList>=obsCentrality])+1)/(nIter+1)
    
    if(returnData){
      
      exRan[[tipo2]] = list(
        obsEpochs = ex2[[tipo2]],
        obsEpochSummary = obsEpochSummary,
        obsCentrality = obsCentrality,
        obsDispersion = obsDispersion,
        
        pvalue = pval,
        epochSummaryFUN = epochSummaryFUN,
        centralityFUN   = centralityFUN,
        dispersionFUN   = dispersionFUN,
        
        cohensList = cohenlist,
        cohensHDR  = hdrcde::hdr(cohenlist,credibilityMass*100)$hdr,
        cohensBest = hdrcde::hdr(cohenlist,credibilityMass*100)$mode,
        
        cliffsList = clifflist,
        cliffsHDR  = hdrcde::hdr(clifflist,credibilityMass*100)$hdr,
        cliffsBest = hdrcde::hdr(clifflist,credibilityMass*100)$mode,
        
        ranCentraList = ranCentraList,
        ranCentraHDR = hdrcde::hdr(ranCentraList,credibilityMass*100)$hdr,
        
        ranDisperList = ranDisperList,
        ranDisperHDR  = hdrcde::hdr(ranDisperList,credibilityMass*100)$hdr,
        
        ranEpochSummary = lapply(resList, \(x){x$randEpochSummary})
      )
      class(exRan) = c("DyadCatPerm", "list")
      attributes(exRan) = c(attributes(exRan), list("parameters" = args))
      
      
    }
    
    


    if(!is.na(plotPrefix)){
      cat("\rInitializing graphical output...")
      
      
      
      #needed functions
      plotHDI <- function( sampleVec , credMass=0.95, y=10, h.len = 5, ...){
        hdi = as.numeric(hdrcde::hdr( sampleVec , credMass*100)$hdr)
        xpar = par()[["xpd"]]
        #there can be multiple intervals in multimodal distributions
        ni = length(hdi)/2
        par(xpd=NA)
        if(ni%%1 != 0){
          #if there is an odd number of hdi intervals, simplify to the broader range
          #but throw a warning
          ni=1
          hdi = c(hdi[1], hdi[length(hdi)])
          warning("Multiple and odd HDR intervals were collapsed to one.")
        }
        for(i in 1:ni){
          x1 = hdi[1+2*(i-1)]; x2 = hdi[2+2*(i-1)]
          segments(x1,y,x2,y,...)
          segments(x1,y-h.len/2,x1,y+h.len/2,...)
          segments(x2,y-h.len/2,x2,y+h.len/2,...)
        }
        par(xpd=xpar)
        
      }
      
      png(paste0(plotPrefix,"_",tipo2,"_",minEpochSec,"s_",nIter,"perm.png"),width = 60*3, height = 60*3, units = "mm",res = 300, type="cairo")
      # par(mfrow=c(1,1),mar=c(5,4,1,1)+0.1,cex=0.9)
      par(cex=0.9, mar=c(5, 4, 4, 4) + 0.1)
      # par(xpd=NA)
      
      if(length(xlim)==1 && xlim=="auto"){
        xmin = min(randDenX, obsCentrality)#, minRandX )
        xmax = max(randDenX, obsCentrality)#, maxRandX)
        
      } else {
        xmin = xlim[1]
        xmax = xlim[2]
      }
      
      breaks = seq(min(min(ranCentraList),xmin),
                   max(max(ranCentraList),xmax),by=(xmax-xmin)/40)
      
      realDen = density(obsEpochSummary,from=xmin,to=xmax)
      maxY = max(realDen$y)
      #histogram of estimated parameters for random
      plotData = hist(ranCentraList, breaks = breaks, plot = F)
      plotData$density = rangeRescale( plotData$density, 0, -maxY)
      toK = which(plotData$density!=0)
      plotData = lapply(plotData, \(x)x[toK])
      # maxY = max(plotData$counts)
      minY = min(plotData$density)
      # maxY = max()
      
      plot(-99999,
           main=paste0(categoryIndex," - ", tipo2,"\nPermutation test on ",nIter," extractions of of ",n, " random epochs each"),
           xlab=paste(epochSummaryFUN, series), ylab = "Density",
           cex.main = 0.9, yaxt="n",
           ylim = c(minY*2,maxY*2),
           xlim=c(xmin,xmax))
      axis(2, at=c(seq(round(minY*2),0) ,seq(0,round(maxY*2))),
           labels = abs(c(seq(round(minY*2),0) ,seq(0,round(maxY*2)))))
      
      polygon(c(xmin,realDen$x,1), c(0,realDen$y,0), col = rgb(0.78, 0.89, 1, alpha = 0.6),border = NA)
      lines(realDen$x, realDen$y, col = rgb(0.58, 0.69, 1, alpha = 1),lwd=2)
      
      polygon(c(xmin,randDenX,1), c(0,-randDenY,0) , col = rgb(0.2, 0.3, 0.2, alpha = 0.4),border = NA)
      drawLines  = if(nIter > nLines) sample(nIter,nLines) else nIter
      cat("\rDrawing lines...                ")
      for(i in 1:length(drawLines)){
        k = drawLines[i]
        lines(resList[[k]]$randDen$x,-resList[[k]]$randDen$y,col=rgb(0.1,0.2,0.1,0.15))
      }
      rect(plotData$breaks[1:40],0,plotData$breaks[2:41],plotData$density,col=0,lwd=3,border = 0)
      rect(plotData$breaks[1:40],0,plotData$breaks[2:41],plotData$density,col=0,lwd=1,border = 1)
      abline(h=0,lwd=2,col=0)
      
      maxY = maxY*1.1
      segments(obsCentrality,0,obsCentrality,maxY,lwd=2,col=2)
      
      text (obsCentrality, maxY+maxY/100*2, paste("p-value =", format(round(pval,4),nsmall = 4) ), cex = 0.8, col=2, font=2 )
      text (obsCentrality, maxY+maxY/100*10, paste0("Observed ", centralityFUN," = ",round(obsCentrality,3)), cex = 0.8, col=2, font=2 )
      
      ## Kruschke, 2014 for .89 HDI
      plotHDI(ranCentraList,lwd=3,y=0, h.len=maxY/100*5, credMass = credibilityMass)
      cohTit = paste0(centralityFUN," effect size:", round(p_centralityFUN(cohenlist),2),
                      " [",credibilityMass,"% Credible Interval:", paste(round(hdrcde::hdr(cohenlist,credibilityMass*100)$hdr,2),collapse = ", "),"]" )
      title(main = cohTit,line = 0.5, cex.main=0.9, font.main=3)
      legend("topleft",
             lty=c(1,0,1,0,0),
             lwd=c(2,1,2,0,0),
             col=c(2,1,1,rgb(0.78, 0.88, 1, alpha = 0.6),rgb(0.4, 0.4, 0.4, alpha = 0.4)),
             legend=c(paste(centralityFUN, "of",n,"real epochs"),
                      bquote("Distribution of the real "*.(centralityFUN)~"under"~H[0]),
                      bquote(.(credibilityMass)*"% Credible Interval of the real "*.(centralityFUN)~"under"~H[0]),
                      "density of observed data",
                      "density of random data"),
             pch = c(NA,0,NA,15,15),pt.cex =c(0,2.5,0,3,3), y.intersp=1.3,
             bty = "n")
      
      
      ###################################################
      dev.off()
    } # end of plot if
    cat("\rDone ;)                         ")   
  } # end of main loop
  
  if(returnData) return(exRan)
  
}



#' @export
print.DyadCatPerm = function (x, ...) {
  pa = attr(x, "parameters")
  cat0("\r\nAfter ",pa$nIter, " permutations of \"",pa$category, " - ",pa$categoryIndex,
       "\" categories,\r\nthe probability of superiority of observed ",pa$centralityFUN,"s to chance,",
       "\r\nunder the null hypothesis, of no temporal association between\r\n",
       pa$series, " and ", pa$category," were:\r\n\r\n")
  c0 = names(x)
  c1 = unlist(lapply(x, \(x)x$pvalue))
  c2 = unlist(lapply(x, \(x)x$cohensBest))
  c3 = unlist(lapply(x, \(x)x$cohensHDR[1,1]))
  c4 = unlist(lapply(x, \(x)x$cohensHDR[1,2]))
  tab = data.frame(c0,c1,c2,c3,c4)
  colnames(tab) = c("categories","p-value","Cohen's d", paste0("d ",pa$credibilityMass*100,"%CI lower"), paste0("d ",pa$credibilityMass*100,"%CI higher"))
  print(tab, row.names=FALSE)
  # 
  # colnames(newline1) = paste0("p_",colnames(newline1))
  # 
  # cat0("\r\nWith the following Cohen's d effect sizes:\r\n\r\n")
  # newline2 = as.data.frame(lapply(x, \(x)x$cohensBest))
  # colnames(newline2) = paste0("ES_",colnames(newline2))
  # print(newline2, row.names=FALSE)
  # 
}
