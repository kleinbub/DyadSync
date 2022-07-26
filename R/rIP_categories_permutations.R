#this file is a revised and controlled (25/07/2022) MASTER file to perform permutation analyses
#specifically, given a well formatted DyadExperiment object with some calculated DyadStream (e.g. synchronization) and DyadCategory data
#it cuts the stream according to the categories intervals, extracts random permutations of the stream, and calculates p-values and effect sizes
#for the observed central values compared to random ones.

#setup
rm(list=ls())
library(DyadSync)
library(hdrcde)

#category settings

#' Permutation tests for signals in specific epochs
#' given a well formatted DyadExperiment object with some calculated DyadStream (e.g. synchronization or physiology) and DyadCategory data
#' it cuts the stream according to the categories intervals, extracts random permutations of the stream, and calculates p-values and effect sizes
#' for the observed central values compared to random ones.
#'
#' @param x a DyadExperiment Object.
#' @param plotPrefix character. path and prefix of output plots names.
#' @param signal character. name of a DyadSignal contained in x.
#' @param sync character. Name of a synchronization object. E.g. PMBest, or CCFBest
#' @param streamKey character. Name of a DyadStream object contained within sync
#' @param category character. Name of a DyadCategory object contained in x
#' @param groupIndex character. Name of a factor column in 'category'.
#' @param absolute logical. should the DyadStream be transformed to absolute values?
#' @param minEpochSec number of seconds of minimum epoch. Shorter epochs will be deleted.
#' @param minOccurrences min occurrences of given category to be kept. Less frequent categories will be deleted.
#' @param keepOnly vector of characters denoting if only given categories must be analyzed. Eg: c("cat1","cat2").
#' @param nIter integer. Number of random extractions.
#' @param rotate if true windows are extracted from a randomly shifted timeseries. This is better TRUE.
#' @param extract_from Either "all", or "remaining". Should new extractions be done randomly from the whole stream or only on non-categorized parts.
#' @param summarizingFunction String. the function used to summarize each epoch's time-series to a single value. Suggested: median.
#' @param parameterFunction  String. the function used to summarize all epochs's single values. Suggested: median.
#'
#' @return
#' @export
#'
#' @examples
categoryPerm = function(x,
                        plotPrefix,
                        signal="SC",
                        sync="PMdev",
                        streamKey="sync",
                        category = "PIRS",
                        groupIndex = "code2",
                        keepOnly = c(), #
                        minEpochSec = 3, #
                        minOccurrences = 10, #
                        
                        absolute = FALSE,
                        
                        #randomization settings
                        nIter = 1000,
                        rotate = TRUE, 
                        extract_from = c("all", "remaining")[2], 
                        summarizingFunction = c("mean","median")[2],
                        parameterFunction = c("mean","median")[2]) {


  d3 = experiment
  if(absolute){
    for(i in 1:length(d3)){
      d3[[i]][[signal]][[sync]][[streamKey]] = abs(d3[[i]][[signal]][[sync]][[streamKey]])
    }
  }
  samplesPerSec = frequency(d3[[i]][[signal]][[sync]][[streamKey]])
  
  d4 = epochStream(d3, signal=signal, sync=sync,streamKey = streamKey,
                   category=category,groupIndex=groupIndex, mergeEpochs = F, artefact.rm=F)
  
  resName = paste0(c(toupper(category), "_", totitle(c(sync,streamKey))),collapse = "")
  ex2 = extractEpochs (d4, signal=signal, epochStream=resName)
  
  
  ## keep only selected categories
  if(length(keepOnly)>0){
    ex2 = ex2[keepOnly]
  }
  ex2 = ex2[sapply(ex2,length)>0]
  cat("\r\nANALYZED CATEGORIES: ", names(ex2))
  
  
  shift <- function(x, n) {
    n = n%%length(x)
    if(n==0) x
    else c(x[(n+1):length(x)], x[1:n])
  }
  
  fastcohen = function(a,b){
    #this requires that a and b don't have missing values
    (mean(a)-mean(b))/sqrt(((length(a)-1)*sd(a)^2 +  (length(b)-1)*sd(b)^2 )/(length(a)+ length(b)-2 ) )
  }
  
  summar = eval(parse(text=summarizingFunction))
  param  = eval(parse(text=parameterFunction))
  
  # na.omit.safe = function
  #remove epochs shorter than minEpochSec seconds
  
  toRemove = c()
  
  for(i in 1:length(ex2)){
    #remove NAs from the sequences, so that they don't influence results
    for(j in 1:length(ex2[[i]])){
      ex2[[i]][[j]] = ex2[[i]][[j]][!is.na(ex2[[i]][[j]])]
    }
    durations = as.numeric(sapply(ex2[[i]],length))
    toKeep = durations > (minEpochSec*samplesPerSec) #e.g. 5 seconds x 10 samplepersecond
    ex2[[i]] = ex2[[i]][toKeep]
    
    #if any category has less than 10 instances, add to remove list
    if(length(ex2[[i]])<minOccurrences) toRemove = c(toRemove, i)
    
  }
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
  for(tipo2 in toPlot){
    # tipo2 = "ASS" ## <-----------------------------------------------####
    print(tipo2)
    #remove NAs
    ex2[[tipo2]] = ex2[[tipo2]][sapply(ex2[[tipo2]], function(x){!all(is.na(x))})] #da FALSE solo se tutta la finestra è NA
    
    #n seq di tipo CODIFICA
    n = length(ex2[[tipo2]])
    
    #1 crea sequenze di finestre consecutive dalle duration
    durations = lapply(ex2, function(x)sapply(x,length))
    durReal = durations[[tipo2]]
    
    #calcola mediana su ciascuna finestra dei dati reali
    epochMediansReal = as.numeric(unlist(lapply(ex2[[tipo2]], summar, na.rm=T)))
    epochMediansReal = epochMediansReal[!is.na(epochMediansReal)]
    
    #calcola statistica centrale sui dati reali
    parReal = param(epochMediansReal,na.rm=T)
    
    if(extract_from =="all") {
      #estrai tutti i segnali sync delle sedute e togli i NA
      full_from = lapply(d4, function(x) as.numeric(x[[signal]][[sync]][[streamKey]]))
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
    
    #setup empy containers
    cohenlist = numeric(nIter)
    parlist = numeric(nIter)
    parlist2 = numeric(nIter)
    res = numeric(sum(durReal))
    madness = matrix(0,nrow=nIter,ncol=512)
    system.time({
      for(k in 1:nIter){
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
        pads = floor(seed/sum(seed)*maxgap)
        if(sum(pads)>maxgap || sum(pads)<maxgap-n) warning("casual paddings had weird length")
        
        #sposta xstart e xend
        rstart = xstart; rend =xend
        for(i in 1:n) {
          rstart[i:n] = rstart[i:n] + pads[i]
          rend[i:n] = rend[i:n] + pads[i]
        }
        rstart
        rend
        sum(rend-rstart)
        
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
          
        }
        
        #estrai il parametro dalle finestre random
        epochMediansRand = as.numeric(unlist(lapply(rfsplit, summar,na.rm=T)))
        epochMediansRand = epochMediansRand[!is.na(epochMediansRand)]
        cohenlist[k] = fastcohen(epochMediansReal,epochMediansRand)
        
        parlist[k]   =  param(epochMediansRand, na.rm=T)
        parlist2[k]  =  sd(epochMediansRand, na.rm=T)
        randDen      =  density(epochMediansRand, from = if(absolute) 0 else -1, to = 1)
        madness[k,]  =  randDen$y
      }
    })
    
    exRan[[tipo2]] = parlist
    
    #needed functions
    plotHDI <- function( sampleVec , credMass=0.95, y=10, h.len = 5, ...){
      hdi = as.numeric(hdr( sampleVec , credMass*100)$hdr)
      x1 = hdi[1]; x2 = hdi[2]
      xpar = par()[["xpd"]]
      par(xpd=NA)
      segments(x1,y,x2,y,...)
      segments(x1,y-h.len/2,x1,y+h.len/2,...)
      segments(x2,y-h.len/2,x2,y+h.len/2,...)
      par(xpd=xpar)
      
    }
    
    png(paste0(plotPrefix,"_",tipo2,"_",minEpochSec,"s_",nIter,"perm.png"),width = 85*2, height = 60*2, units = "mm",res = 300, type="cairo")
    # par(mfrow=c(1,1),mar=c(5,4,1,1)+0.1,cex=0.9)
    par(cex=0.9, mar=c(5, 4, 4, 4) + 0.1)
    par(xpd=NA)
    
    ## big nice plot
    breaks = seq(-1,1,by=0.05)
    plotData = hist(parlist, breaks = breaks, plot = F)
    maxY = max(plotData$counts)
    # maxY = max(plotData$density)
    
    xmin = min(plotData$breaks,parReal)
    xmax = max(plotData$breaks,parReal)
    if(absolute) xmin = 0 else xmin = -1
    
    hist(parlist,breaks = breaks,
         main=paste0(groupIndex,"-", tipo2,"\nDistribution of ",nIter," ",parameterFunction,"s of ",n, " random epochs each"),
         xlab=paste(summarizingFunction, "synchronization"),
         cex.main = 0.9, col=rgb(1,1,1,1),
         ylim = c(0,maxY*1.2),
         # xlim = c(xmin,xmax))
         xlim=c(xmin,1))
    
    realDen = density(epochMediansReal,from=if(absolute) 0 else -1,to=1)
    polygon(c(xmin,realDen$x,1),c(0,realDen$y,0)*maxY, col = rgb(0.78, 0.89, 1, alpha = 0.6),border = NA)
    polygon(c(xmin,randDen$x,1), c(0,apply(madness,2,median)*maxY,0) , col = rgb(0.4, 0.4, 0.4, alpha = 0.4),border = NA)
    hist(parlist,breaks = breaks,add=T, col=rgb(1,1,1,0.6))
    
    pval =(length(parlist[parlist>=parReal])+1)/(nIter+1)
    segments(parReal,0,parReal,maxY,lwd=2,col=2)
    # segments(mean(parlist),0,mean(parlist),maxY/2,lwd=2,col=rgb(0.4,0.4,0.4))
    
    text (parReal, maxY+maxY/100*2, paste("p-value =", format(round(pval,4),nsmall = 4) ), cex = 0.8, col=2, font=2 )
    text (parReal, maxY+maxY/100*6, paste0("Observed ", parameterFunction," = ",round(parReal,3)), cex = 0.8, col=2, font=2 )
    
    axis(4, at = seq(0, maxY*1.2, length.out=7), labels = round(rangeRescale(seq(0, maxY*1.2, length.out=7), 0, max(realDen$y)),1))
    mtext("Density", side = 4,line=2.8)
    
    
    # text(parReal,510, "real data",col=2)
    # text(median(parlist),-2,"95% HDI",cex=0.9,font=2)
    
    plotHDI(parlist,lwd=3,y=0, h.len=5)
    cohTit = paste0(parameterFunction," effect size:", round(param(cohenlist,na.rm=T),2),
                    " [89% HDR range:", paste(round(hdr(cohenlist,89)$hdr,2),collapse = ", "),"]" )
    title(main = cohTit,line = 0, cex.main=0.9, font.main=3)
    legend("topleft",lty=c(1,0,0),lwd=c(2,0,0),col=c(2,rgb(0.78, 0.89, 1, alpha = 0.6),rgb(0.4, 0.4, 0.4, alpha = 0.4)),
           legend=c(paste(parameterFunction, "of",n,"real epochs"), "density of observed data", "density of random data"),
           pch = c(NA,15,15),pt.cex =c(0,3,3), y.intersp=1.3,
           bty = "n")
    
    
    ###################################################
    dev.off()
    
  }




}