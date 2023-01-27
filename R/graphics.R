###NOT @ratsOK ç_ç

##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  
## DyadPlot.R
# Main TS plotting tool
#
############################################################################################################ 
## Changelog
# v1.3 - the world is burning. New plots for ppBest/ppSync. Everything must be checked.
# v1.1 - added boxes plotting
# v1.0 stable
#
############################################################################################################ 
## ToDo
# -IMPORTANT: CCF & best lag must NOT be centered!! (works only with scale="normalize")
# -legends
# -automatic csv reader as boxes
# -color code connection lines
############################################################################################################ 
## Credits
# Author: Johann R. Kleinbub
# Contact: johann.kleinbub@gmail.com
############################################################################################################ 

# NB the traditional colors of the lines are:
## deeppink3
## dodgerblue3


##Main TS plotting tool

#' Easy axes for ts objects
#'
#' @param x a ts object
#' @param side side of the axie
#' @param every interval between ticks, in seconds
#' @param out format of the time, one of "hour","min","sec" for respectively the formats
#' hh:mm:ss, mm:ss, ss.
#' @param las 'label axis style' sets the labels orientation. See ?par for details
#' @param ... further arguments passed to axis
#'
#' @return
#' @export
#'
#' @examples
tsAxis = function(x, side, every=1, out=c("hour", "min", "sec"), las=1,...){
  out = match.arg(out,c("hour", "min", "sec") )
  axis(side, at=as.numeric(time(x)[c(T,rep(F,every*frequency(x)-1))]), labels = timeMaster(round(time(x)[c(T,rep(F,every*frequency(x)-1))]),out),las=las, ... )
  
}

#' @export
#'
plot.DyadExperiment =function(x, signal, ...){
  print("plot.DyadExperiment")
  if(missing(signal)){
    res = unique(unlist(lapply(x,names)))
    if (length(res)>1) stop ("multiple signals detected, please specify which you want to print")
    else signal = res
  }
  for(i in seq_along(x)){
    plot(x[[i]][[signal]])
  }
}

#' @export
#'
basicPlot = function(x, signal, ...){
  UseMethod("basicPlot",x)
}
#' @export
#'
basicPlot.DyadExperiment =function(x, signal, ...){
  print("plot.DyadExperiment")
  if(missing(signal)){
    res = unique(unlist(lapply(x,names)))
    if (length(res)>1) stop ("multiple signals detected, please specify which you want to print")
    else signal = res
  }
  for(i in seq_along(x)){
    basicPlot(x[[i]][[signal]], main=paste(dyadId(x[[i]]),"-",signal ))
  }
}
#' @export
#'
basicPlot.DyadSignal = function(x, ...) {
  print("plot.DyadSignal")
  par(mfrow=c(2,1), mar=c(0,2.5,2.5,0), cex = 0.8)
  plot((x$s1),xaxt="n" ,ylab="uS", ...)
  par(mar=c(3.5,2.5,0,0))
  plot(time(x$s2),x$s2,col="red",xaxt="n",t="l",xlab="", cex.main = 0.001, ...)
  tSteps = round(time(x$s2)[seq(1, length(time(x$s2)),by=frequency(x)*60 )])
  axis(1,at = tSteps, labels = timeMaster(tSteps,"min"),tick = T,las=2 )
}






#generate as different as possible colors
#' Title
#'
#' @param n 
#' @param demo 
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
mycolz = function(n,demo=F,alpha=1){
  fmod= function(k,m){
    j = floor(k/m)
    a = k-m*j
    return(a)
  }
  colz =numeric(n)
  for (i in 1:n) {
    colz[i] = hsv(fmod(i * 0.618033988749895, 1.0),
                  0.8, sqrt(1.0 - fmod(i * 0.618033988749895, 0.5)),alpha = alpha)
    
  }
  if(demo)pie(rep(1,n), col=colz)
  return(colz)
}
dual.colors=function (n, alpha = 1, center=0, neg.hue=0, pos.hue=0.4, stretch = 1, demo=F)
  #good hue values
  #min 0.05 max 0.6
  #min 0.05 max 0.4
  #stretch = 1 indica una distribuzione lineare dei colori, da 2 in su è esponenziale
  #v1.1 - aggiunto "center" che permette di spostare il bianco dallo zero ad un centro
  #numericamente corretto
  #     - aggiunto un intervallo attorno a mean.hue per evitare la zona "gialla uguale in alto e in basso" 
{
  if ((n <- as.integer(n[1L])) > 0L) {
    if(center<1 && center>-1) center = center*100
    even.n <- n%%2L == 0L
    k <- n%/%2L
    l1 <- k + 1L - even.n +center #sequenza di numeri minori di center
    l2 <- n - k + even.n -center #sequenza di numeri maggiori di center
    mean.hue = mean(c(neg.hue, pos.hue))
    res = c(
      if (l1 > 0L)
        hsv(
          h = rangeRescale( seq(mean.hue*0.65,neg.hue, length.out = l1) ^stretch , mean.hue*0.65,neg.hue),
          s = seq.int(0.7, ifelse(even.n,0.5/k, 0), length.out = l1),
          v = seq(0.8, 0.9, length.out = l1), alpha = alpha), 
      if (l2 > 1)
        hsv(
          h = rangeRescale( seq(mean.hue*1.35,pos.hue,length.out = l2)[-1L] ^stretch , mean.hue*1.35,pos.hue),
          s = seq.int(0, 0.7, length.out = l2)[-1L], 
          v = seq.int(0.9, 0.8, length.out = l2)[-1L], alpha = alpha)
    )
  }
  else res = character()
  if(demo){
    plot(1,type="n",xlim=c(0,4),ylim=c(-100,100))
    boxh=(200/(n-1))/2
    rect(1,   seq(-100,100,length.out = n)-boxh,  3,  seq(-100,100,length.out = n)+boxh, col=res,border=NA)
  }
  res
}


#rip of cm.colors()
dual.colors.deleteme=function (n, alpha = 1, neg.hue=0, pos.hue=0.4,demo=F)
  #good hue values
  #min 0.05 max 0.6
  #min 0.05 max 0.4
{
  if ((n <- as.integer(n[1L])) > 0L) {
    even.n <- n%%2L == 0L
    k <- n%/%2L
    l1 <- k + 1L - even.n
    l2 <- n - k + even.n
    mean.hue = mean(c(neg.hue, pos.hue))
    res = c(
      if (l1 > 0L)
        hsv(
          h = seq(neg.hue, mean.hue, length.out = l1),
          s = seq.int(0.7, ifelse(even.n,0.5/k, 0), length.out = l1),
          v = seq(0.8, 1, length.out = l1), alpha = alpha), 
      if (l2 > 1)
        hsv(
          h = seq(mean.hue,pos.hue,length.out = l2)[-1L],
          s = seq.int(0, 0.7, length.out = l2)[-1L], 
          v = seq.int(1, 0.8, length.out = l2)[-1L], alpha = alpha)
    )
  }
  else res = character()
  if(demo){
    plot(1,type="n",xlim=c(0,4),ylim=c(-100,100))
    boxh=(200/(n-1))/2
    rect(1,   seq(-100,100,length.out = n)-boxh,  3,  seq(-100,100,length.out = n)+boxh, col=res,border=NA)
  }
  res
}
# poi=dual.colors(101,demo=T)
# poi=dual.colors(101,demo=T,neg.hue = 0,pos.hue = 0.4)

catBoxes = function(categ, column, cols=mycolz(length(levels(categ[[column]])),alpha=0.5) ){
  df = categ[c("start","end",column)]
  df$color = cols[as.integer(df[[column]])]
  return(df)
}


#' Plots a window of signal and the AMICo algorithm results
#'
#' @param signal object of class DyadSignal
#' @param start start of plotting window. "mm:ss" or a number of seconds
#' @param end end of plotting window. "mm:ss" or a number of seconds
#'
#'
plotAMICo = function(signal, start, end, syncName = "AMICo"){
  #if()
  rs1x = rangeRescale(signal$s1,0,1.5)-0.6
  rs2x = rangeRescale(signal$s2,0,1.2)-0.25
  
  xStart = timeMaster(start,"s")
  xEnd   = timeMaster(end,"s")
  
  rs1 =  window(rs1x,start=xStart,end=xEnd)
  rs2 =  window(rs2x,start=xStart,end=xEnd)
  k1 = min(rs1-mean(rs1))+mean(rs1)
  k2 = min(rs2-mean(rs2))+mean(rs2)
  rs1 = rs1-k1
  rs2 = rs2-k2
  
  time = window(signal$time, start=xStart,end=xEnd)
  valid = ts(rep(T, length(res1) ), frequency = frequency(rs1),start=start(rs1))
  for(i in 1:nrow(signal$artefacts)){
    window(valid, start=signal$artefacts$start[i],end=signal$artefacts$end[i]) <- F
  }
  
  sync1 = rangeRescale(window(signal[[syncName]]$sync,start= xStart, end=xEnd), 0.6,0.8,-1,1)

  plot(rs1, col=0, lty=3,xaxs="i",xaxt="n",yaxt="n",ylim=c(-0.05,0.81),xlab="",ylab="")
  lines(time[valid], rs1[valid],col="#00A479",lwd=2)
  lines(time[valid], rs2[valid],col="#3173BA",lwd=2,lty=1)
  xbest2 = signal$AMICo$xBest[signal[[syncName]]$xBest$syncBound==T,]
  segments(time(signal$s1)[xbest2$s1],rs1x[xbest2$s1]-k1,
           time(signal$s2)[xbest2$s2],rs2x[xbest2$s2]-k2,col=1,lwd=2)
  
  lines(sync1,lty=1,col="#DCAA27",lwd=2) 

  abline(h=c(0.6,0.8,0.70),lty=c(1,1,3))
  tSteps = signal$time[seq(1, length(signal$time),by=100 )]
  axis(1,at = tSteps, labels = timeMaster(tSteps+1,"min"),tick = T,las=1 )
  axis(2, at=c(0.6,0.8,0.70),labels = c(-1,1,0),las=2)
  mtext("Correlation",2,las=3,at = 0.70,line = 2,col=1)
  mtext("Normalized\nskin conductance",2,las=3,at = 0.200,line = 0.5,col=1)
  title(xlab="Time (mm:ss)")
  
}



plotSignal = function (lineList, path="test.svg", boxList=NULL, connect=T, stacked=F, videoSec=0,
                       preview=T, preview_time=NULL, preview_length=60,
                       scale = c("corRange","normalize","raw","windowed"),...){

  dots <- list(...) 
  ll = lineList
  #fixed settings. These are just good
  #w=600 #lol invece fallo dinamico sulla durata
  h=15 
  
  halfBoxH = 0.5 #height of boxes
  winScaleSec = 15 #for windowed scale, cut every 15 secs.
  scale = match.arg(scale) 
  # mynames = colnames(dyad)
  # xyl = list( dyad[[1]], #patient
  #             dyad[[2]]  #clinician
  # )
  
  #are all frequencies equal?
  if(!diff(range(sapply(ll,frequency))) < .Machine$double.eps ^ 0.5){
    print("signal have different frequencies, so what?")
  }
 


  
  # if(preview){
  #   #select 2 random minutes of signal
  #   min_v = min(sapply(ll,length))*0.1
  #   max_v = min(sapply(ll,length))*0.9
  #   #ran_start= round(runif(1,min_v,max_v))
  #   ran_start=min_v
  #   #ll = lapply(ll, function(x){ts(x[1:60*frequency(x)+ran_start],frequency = 1, start = ran_start/(60*frequency(x)) )})
  # }
  if(preview){
    min_v = ceiling(min(sapply(ll, end)[1,])*0.1)
    max_v = ceiling(min(sapply(ll, end)[1,])*0.8)
    if(is.null(preview_time))
      ran_start= round(runif(1,min_v,max_v))
    else ran_start = preview_time
    #ran_start=min_v
    ll = lapply(ll, function(x){ window(x, start=ran_start, end=ran_start+preview_length)})
    
  }
  
  ##normalize the signal
  if(scale=="normalize") {
    ll = lapply(ll, function(x){
      if(grepl("CCF",attr(x,"name"))){
        x
      } else {
        scale(x)
      }
    })
  } else if (scale=="windowed") {ll = lapply(ll, stepCenter, winSec=winScaleSec)
  } else if (scale=="corRange") {ll = lapply(ll, rangeRescale, rangeMin=-1, rangeMax =1)}
  #   sync = ts(syncMat[,(ncol(syncMat)-1)/2+1], frequency= 1/winSec) #default è lag 0
  

  begin = min(sapply(ll, start)[1,]) #earliest observation, in seconds
  duration = max(sapply(ll, end)[1,]) #max duration, in seconds
  w = duration* 0.15
  #if(is.null(lineCol)) colz = mycolz(10,F) else colz = lineColz
  if(!preview)  svg(path, width = w, height= h)
  if(!"ylim" %in% names(dots)) myYlim = c(min(sapply(ll,min,na.rm=T),-1.2),max(sapply(ll,max, na.rm=T),1.2)+ifelse(stacked,5,0))
     else myYlim = dots$ylim
  
  # par(mar=c(5,0,4,0)+0.1)
  plot(ll[[1]],lwd=2,xaxp=c(0,duration, duration),xlim=c(begin,duration), ylim=myYlim, type="n", xaxt='n', bty="n",xaxs = "i")
  
  ##crea l'asse x coi minuti
  tickSec=round(begin):round(duration) + videoSec
  tickMin = floor(tickSec/60)
  secLeft = tickSec-(tickMin*60)
  #aggiungi gli zeri
  tickMin[tickMin<10] = paste0('0',tickMin[tickMin<10])
  secLeft[secLeft <10] = paste0('0', secLeft[secLeft<10])
  ticks = paste(tickMin,secLeft, sep=":")
  if(preview){
    axis(1, at=(round(begin):round(duration))[c(T,F)], labels=ticks[c(T,F)],las=2,pos=myYlim[1]-0.2 )
  } else axis(1, at=round(begin):round(duration), labels=ticks,las=2,pos=myYlim[1]-0.2 )
  
  
  yModifier = rep(0,length(ll))
  if(stacked) yModifier[clinician] = 5 
  
  #####TESTAREA
  if(connect){
    #downsample high sampling rate signals
    ll2 = lapply(ll, function(s){
      if(frequency(s)>5) {
        cat("\r\nDecimating ",attr(s,"name"))
        streamDecimate(s,newSampRate = 5)
      } else s
    })
    minLen = min(sapply(ll2, length))
    ll2 = lapply(ll2, function(x){
      window(x, start=start(x), end=c(trunc(minLen/frequency(x)), (minLen/frequency(x))-trunc(minLen/frequency(x)) ))
      })
    
    patient =   which(grepl("patient",   sapply(ll,attr,"name"), ignore.case = T))[1]
    clinician = which(grepl("clinician", sapply(ll,attr,"name"), ignore.case = T))[1]
    if(is.na(patient) || is.na(clinician)) stop("patient or clinician data could not be identified, please rename signals or set 'connect' to FALSE")
    ccfval = which(grepl("bestCCF|lag-*\\d|ppCor", sapply(ll,attr,"name"), ignore.case = T))[1]
    if(is.na(ccfval)) stop("no signal with name attribute == 'bestCCF' or lagXX.")
    if(grepl("lag(-*\\d+)", name(ll[[ccfval]]))) {
      blag = as.numeric(gsub("(.*lag)(-*\\d+)(.*)","\\2", name(ll[[ccfval]])))
      ll2[[length(ll2)+1]] = rep(blag, length(ll2[[ccfval]]))
      blag = length(ll2)
    } else {
      blag = which(grepl("bestLag|ppLag", sapply(ll,attr,"name"), ignore.case = T))[1]
      
    } 
    if(is.na(blag))stop("no signal with name attribute == 'bestLag' and no constant lag found.")
      
      
    
    
    y1_coords = (time(ll2[[clinician]]) + ll2[[blag]]-begin)*frequency(ll[[clinician]])
    
    y1_coords[y1_coords<1 ]= NA
    
    all_col= dual.colors(201)
    iCCF = as.vector(round(ll2[[ccfval]],2)*100 + 101)
    


    
    # plot(ll[[patient]]) #paziente
    # lines(ll[[clinician]],col="blue") #terapeuta
    segments(x0= time(ll2[[patient]]), y0 = ll2[[patient]],
             x1=time(ll2[[clinician]]) + ll2[[blag]],
             y1= ll[[clinician]][y1_coords]+yModifier[clinician] ,col=all_col[iCCF]
    )
    
    #disegna i segnali [tranne blag e ccfval]
    mapply(function(lin, i){
      print(i)
      lines(lin+yModifier[i], col=attr(lin,"col"),lwd=attr(lin,"lwd"), lty=attr(lin,"lty"))
      }, ll[-c(blag,ccfval)], seq_along(ll)[-c(blag,ccfval)],SIMPLIFY=T)
    
    
  } else {
    #disegna i segnali
    
    #questo li stacka uno sopra l'altro ma se li voglio sovrapposti?
    # mapply(function(lin, i){lines(lin+(i-1)*5, col=attr(lin,"col"),lwd=attr(lin,"lwd"), lty=attr(lin,"lty"))}, ll, seq_along(ll),SIMPLIFY=T)
    mapply(function(lin, i){
      lines(lin+yModifier[i], col=attr(lin,"col"),lwd=attr(lin,"lwd"), lty=attr(lin,"lty"))
      }, ll, seq_along(ll),SIMPLIFY=T)
    
  }
 
  
  #disegna le scatole
  halfBoxH = abs(myYlim[2] - myYlim[1])*0.4
  for(boxes in  boxList)
  apply(boxes,1,function(x){
    rStart = as.numeric(x["start"])
    rEnd   = as.numeric(x["end"]) 
    rect(rStart, 0-halfBoxH,rEnd,0+halfBoxH, col=x["color"], border="black",lty=2)
    text(rStart+(rEnd - rStart)/2-nchar(x[3])*0.1, 0+halfBoxH+0.1, labels=as.character(x[3]),cex=2 )
  })
  
  ccf = attributes(ll[[which(grepl("bestCCF", sapply(ll,attr,"name"), ignore.case = T))[1]]])
  title("SyncPlot")
  mtext(paste0("lagSec: ",ccf$lagSec," incSec: ",ccf$incSec," winSec: ",ccf$winSec," accelSec: ",ccf$accelSec," weight: ",ccf$weight," interpolated: ",ccf$interpolated),line=0.2)
  

  # plot(ll[[1]])
  # lines(ll[[2]],col="blue")
  # segments(x0= time(ll[[2]]), y0 = ll[[2]],
  #   x1=time(ll[[1]])+ll[[7]],
  #   
  #   y1=ll[[1]][time(ll[[2]])+ll[[7]]]
  # 
  #   ,col="blue"
  # )
  # 
  ## nice!!
  # plot(ll[[1]])
  # lines(ll[[2]],col="blue")
  # segments(x0= time(ll[[1]]), y0 = ll[[1]],
  #   x1=time(ll[[2]])+ll[[7]],
  #   #y1=ll[[2]][time(ll[[2]])+ll[[7]]]
  #   y1 = ll[[1]]
  # )

  
  ####

  #linee orizzontali
  
  abline(h=c(-2,-1,0,1,2), lty=3)
  # if(!is.null(boxes[[1]]) & !is.null(boxes[[2]]) ){
  #   if(length(boxes[[1]]) != length(boxes[[2]])) stop("Boxes start and end must have same length")
  #   ##draw boxes!
  #   if (is.null(boxes[[3]])) boxes[[3]] = rep(rgb(0.5,0.5,0.5,0.2), length(boxes[[1]]))
  #   if (is.null(boxes[[5]])) boxes[[5]] = rep(0, length(boxes[[1]]))
  #   
  #   for(i in 1:length(boxes[[1]])){
  #     rect(boxes[[1]][i] +1,boxes[[5]][i]-halfBoxH,boxes[[2]][i]+1,boxes[[5]][i]+halfBoxH, col =boxes[[3]][i], border="black",lty=2)
  #     text((boxes[[2]][i]-boxes[[1]][i])/2+boxes[[1]][i],boxes[[5]][i]+halfBoxH+0.1, labels=as.character(boxes[[4]][i]),cex=2 )
  #   }
  # }
  if(scale=="windowed") {
    abline(v=seq(begin,duration, by=winScaleSec),lty=3,col="darkgrey")
  }
  # syncNames = Map(function(sy,win,inc){paste("Cross-cor with",win,"sec. windows and",inc,"sec. increments.")}, syncList, winSecList,incSecList)
  # print("syncNames:")
  # print(syncNames)
  # legend("topleft", legend=c(mynames,syncNames), col=c(colz[1:2],heat.colors(length(syncList))),lty=c(1,1, rep(3,length(syncList))), lwd=2, cex=2)
  
  box()
  
  if(!preview) graphics.off()
  
}

ppBestPlot = function (signal,savePath){
  #sta roba dovrai integrarla in plotSignal con tutti i crismi
  if (TRUE){
    w = end(signal$s1)[1]* 0.15
    #if(is.null(lineCol)) colz = mycolz(10,F) else colz = lineColz
    svg(savePath, width = w, height= 15)
    sampRate = signal$sampRate
    xbest = signal$ccf$ppBest
    # startx = round(runif(1,1,trunc(length(signal$s1)/signal$sampRate )-preview_sec))
    # xbest = xbest[xbest$s1>startx*signal$sampRate & xbest$s1<(startx+preview_sec)*signal$sampRate,]
    
    # d = window(signal$s2,start=startx, end=(startx+preview_sec)) #[primi 2 min di segnale]
    # d2 = window(signal$s1,start=startx,  end=(startx+preview_sec)) #[primi 2 min di segnale]
    
    d =signal$s2
    d2 =signal$s1
    videoSec = start(d)[1]

    
    sgol_p=2  #l'ordine dei polinomi del filtro
    sgol_n=25 # il numero di sample in ogni finestra: 55 trova i macropicchi, 35 trova tutti i picchi
    fd5  = ts(sgolayfilt(d,  p =sgol_p, n = sgol_n),frequency = 10,start=start(d))
    fd52 = ts(sgolayfilt(d2, p =sgol_p, n = sgol_n),frequency = 10,start=start(d))

    
    par(mfrow=c(1,1))
    kd = rangeRescale(fd5,0,1)
    kd2 = rangeRescale(fd52,0,1)
    plot(kd, ylim=c(0,2),type="n", xaxt='n', bty="n",xaxs = "i",
         main=paste("[PPA v1.93] lag:",signal$ccf$settings$lagSec),lwd=2,col=attr(signal$s2,"col"))
    # title(sub = "z-normalized signal",line=-27.5)
    
    
    begin = min(start(signal$s1)[1]) #earliest observation, in seconds
    duration = max(end(signal$s1)[1]) #max duration, in seconds


    ##crea l'asse x coi minuti
    tickSec=round(begin):round(duration) + videoSec
    tickMin = floor(tickSec/60)
    secLeft = tickSec-(tickMin*60)
    #aggiungi gli zeri
    tickMin[tickMin<10] = paste0('0',tickMin[tickMin<10])
    secLeft[secLeft <10] = paste0('0', secLeft[secLeft<10])
    ticks = paste(tickMin,secLeft, sep=":")
    axis(1, at=round(begin):round(duration), labels=ticks,las=2,pos=+0.2 )
    
    ##construction grids on s1
    # abline(v=time(d)[xbest$s1])
    # lines(x = seq(time(d)[xbest$s1[1]] , time(d)[xbest$s1[nrow(xbest)]],by=0.1),
    #       y = rangeRescale(c(rep(xbest$sync[1:(nrow(xbest)-1)], diff(xbest$s1)),NA),0,1,-1,1),
    #       col=2)
    # abline(h=c(1,0,0.5), lty=3)
    # dev.off()

    #da qui la parte  di CONNECT
    downsamp = round(seq(1, length(d),by=sampRate/2)) #per chiarezza plotta solo ogni mezzo secondo
    s1Time = time(d)[downsamp] #coordinate temporali di s1 (le x sul grafico blu)
    s2Time = time(d2)[downsamp] + signal$ccf$ppSync$lag[downsamp]/sampRate #per s2 aggiungi il lag (trasformato in secondi)
    s2y = downsamp + signal$ccf$ppSync$lag[downsamp]
    s2y[s2y<1]=NA
    
    all_col= dual.colors(100,neg.hue = 0, pos.hue = 0.4, stretch = 5)
    # all_col= c(heat.colors(100),0, rev(heat.colors(100) ))
    iDD = as.vector(round(signal$ccf$ppSync$sync[downsamp],2)*100 )
    iDD[iDD<0] = 1
    segments(s1Time, kd[downsamp],
             s2Time, kd2[s2y]+1,
             col = all_col[iDD]
    )
    #try with flat shades
    all_col= dual.colors(201,neg.hue = 0, pos.hue = 0.4, stretch = 5)
    for (i in 1:(nrow(xbest)-1) ){
      polygon( x = c(time(d)[xbest$s1[i]:xbest$s1[i+1]], rev(time(d)[xbest$s2[i]:xbest$s2[i+1]])),
               y = c(kd[xbest$s1[i]:xbest$s1[i+1]], rev(kd2[xbest$s2[i]:xbest$s2[i+1]]+1)),
               col = all_col[round(signal$ccf$ppBest$sync[i],2)*100+101], border=NA
      )
    }
    
    #linee nere di connection picco-picco 
    for (i in 1:nrow(xbest)){
      segments(time(d)[xbest$s1[i]],kd[xbest$s1[i]],
               time(d)[xbest$s2[i]],kd2[xbest$s2[i]]+1
      )
    }
    lines(kd2+1, col=attr(signal$s1,  "col"),lwd=2)
    lines(kd   , col=attr(signal$s2,"col"),lwd=2)
    
    

    dev.off()
  }
}

rescaleByWin = function(x, winSec, rangeMin, rangeMax){
  win = winSec*frequency(x)
  win2 = trunc(win/2)
  len = length(x)
  x2 = x
  for(i in seq_along(x2)){
    
    a = max(1,i-win2)
    b = min(i+win2, len)
    if(i<=win2) pos = i else pos = win2
    
    x2[i] = rangeRescale(x[a:b], rangeMin, rangeMax)[pos]
  }
  x2
}



#' Plot a DyadSignal object
#'
#' @param x a DyadSignal object
#' @param sync Either NA or a string pointing to the name of a PMBest or CCFBest object
#' @param ... 
#'
#' @export
#'
plot.DyadSignal = function(x, sync=NA, rescale = c("none","win", "fixed"), ...) {
  if(missing(sync)) sync = NA
  if(missing(rescale)) rescale = "fixed"
  
  dots = list(...)
  rescale = match.arg(rescale, c("none","win", "fixed") )
  if(rescale=="none"){
    rs1=x$s1
    rs2=x$s2
    myYlim = c(min(c(rs1,rs2),na.rm = T), max(c(rs1,rs2),na.rm = T))
    
  } else if(rescale == "win"){
    rs1= rescaleByWin(x$s1, winSec = 60*5,rangeMin = 0,rangeMax = 1) #experimental amplification
    rs2= rescaleByWin(x$s2, winSec = 60*5,rangeMin = 0,rangeMax = 1)
    myYlim = c(0,1)
  } else if(rescale == "fixed"){
    rs1 = rangeRescale(x$s1,0,1, quantile(x$s1,0.025),quantile(x$s1,0.975))
    rs2 = rangeRescale(x$s2,0,1, quantile(x$s2,0.025),quantile(x$s2,0.975))
    myYlim = c(0,1)
    
  }else {
    
  } 
  
  if(!missing(sync) && !is.na(sync)){
    myYlim[2] = myYlim[2] * 1.7
  } else  myYlim[2] = myYlim[2] * 1.1

  dots = setArg("ylim",myYlim, dots)
  dots = setArg("lty", 3, dots)
  dots = setArg("xaxs","i", dots)
  dots = setArg("xaxt","n", dots)
  dots = setArg("ylab",name(x), dots)


  do.call("plot.ts",c(list("x"=rs1),dots))
  # plot(rs1,t="l")
  lines(rs2,lty=3)
  x$time = as.numeric(time(x$s1))
  
  #colorize the good x
  lines(x$time, rs1,col=attr(x$s1,"col"))
  lines(x$time, rs2,col=attr(x$s2,"col"))
  
  tSteps = round(x$time[seq(1, length(x$time),by=sampRate(x) )])
  axis(1,at = tSteps, labels = timeMaster(tSteps,"min"),tick = T,las=2 )
  
  if(!missing(sync) && !is.na(sync)){
    if(class(x[[sync]]) == "PMBest"){
      #draw peak matching connections
      segments(time(x$s1)[x$PMdev$xBest$s1],rs1[x$PMdev$xBest$s1],
               time(x$s2)[x$PMdev$xBest$s2],rs2[x$PMdev$xBest$s2],lty=3)
      xbest2 = x$PMdev$xBest[x$PMdev$xBest$syncBound==T,]
      segments(time(x$s1)[xbest2$s1],rs1[xbest2$s1],
               time(x$s2)[xbest2$s2],rs2[xbest2$s2])
    }

    
    sync = x$PMdev$xBest[x$PMdev$xBest$syncEnd !=0,]
    sync$sync = rangeRescale(sync$sync,0,2,-1,1)
    sync$sync = sync$sync^2
    sync$sync = rangeRescale(sync$sync,1.1,1.6,0,4)
    
    colors= c("#f70c3f","#c16107", "#878787","#878787","#878787", "#86a817", "#23c647")
    colfunc <- grDevices::colorRampPalette(colors=colors, bias=1)
    legendSteps =20
    colz = colfunc(legendSteps)
    colz = c(colz,"#aaaaaa") #colore degli NA
    
    bins = seq(0.5/legendSteps,0.5, by=0.5/legendSteps)+1.1
    
    iCol = sapply(sync$sync,function(x) sum(x>bins)+1)
    iFill = iCol
    iCol[is.na(iCol)]=length(colz)
    iFill[!is.na(iFill)] = NA
    iFill[is.na(iFill)] = 40
    abline(h=quantile(sync$sync))
    abline(h=median(sync$sync),lwd=2)
    rect(round(time(x$s1)[sync$syncStart]),1.1,round(time(x$s1)[sync$syncEnd]),sync$sync,col=colz[iCol])
    
    text(start(rs1)[1]+3,y=median(sync$sync)+0.02,labels = "median sync")
    text(start(rs1)[1]+3,y=quantile(sync$sync,0.25)+0.02,labels = "25%")
    text(start(rs1)[1]+3,y=quantile(sync$sync,0.75)+0.02,labels = "75%")
    
    abline(h=quantile(sync$sync),lty=3)
  }
  
}

#' Plot characters
#' 
#' if it contains a colors list :-)
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.character = function(x, ...){
  if(all(substr(x,1,1)=="#")){
    #these are colors!
    n = length(x)
    nr =1 
    nc= n
    if(n > 5) {
      nr=ceiling(n/5)
      nc =5
    }
    mata = matrix(1:(nc*nr),ncol=nc, byrow=T)
    oldPar = par()
    par(mar=c(0,0,0,0))
    plot(-1000, xlim=c(0,nc),ylim=c(nr,0), ...)
    rect(xleft   = 1:nc-1,
         xright  = 1:nc,
         ytop    = rep(1:nr-1, each=nc),
         ybottom = rep(1:nr  , each=nc),
         col=x, border="#ffffff"
         )
    text( x = 1:nc-0.4,
          y = rep(1:nr-0.85, each=nc),
          labels = c(t(mata)),
          col="#ffffff"
    )
    text( x = 1:nc-0.6,
          y = rep(1:nr-0.85, each=nc),
          labels = c(t(mata)),
          col="#000000"
    )
    text( x = 1:nc-0.5,
          y = rep(1:nr-0.25, each=nc),
          labels = x,
          col="#ffffff"
    )
    text( x = 1:nc-0.5,
          y = rep(1:nr-0.15, each=nc),
          labels = x,
          col="#000000"
    )
    par = oldPar
    
  } else {
    NextMethod()
  }
}
