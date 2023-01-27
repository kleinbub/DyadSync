#' v2RC0.4
#' This is just a packetized version of Amico_2_v4.R
#' please go there to make edits
#' 

# lr,signal = "SC",lagSec=4,match_threshold=0.25,
# minSizeSec=5,weightMalus = 35,algorithm = "v1.1",outputName = "PMdev",minPeakDelta = 0.05

#' Title
#'
#' @param signal 
#' @param lagSec 
#' @param weightMalus 
#' @param match_threshold 
#' @param minSizeSec 
#' @param outputName 
#' @param interval_sec 
#' @param maxPeakDuration 
#' @param ... Further arguments passed to peakFinder
#'
#' @return
#' @export
#'
#' @examples
AMICo2 = function(signal, lagSec, 
                  weightMalus = 50, match_threshold = 0.0, minSizeSec=4,
                  outputName = "AMICo2", maxPeakDuration = "hrt",
                  interval_sec = 15, sgol_p,sgol_n,correctionRangeSeconds,minPeakAmplitude )
{
  args <- c(as.list(environment()), #list(...),
            list(
              sessionId=sessionId(signal),
              dyadId=dyadId(signal),
              groupId=groupId(signal))
  )
  args["signal"] = NULL
  

  
  
  # stop("debug")
  # lagSec=4;
  # sgol_p = 2; sgol_n = 25;  weightMalus = 35;
  # match_threshold = 0.25; minSizeSec=1;
  # algorithm=c("v1.1");
  # outputName = "PMBest_test";
  # correctionRangeSeconds = 0.5; minPeakAmplitude = 0.05;
  # iSession = 1
  # debug = TRUE
  # cols1 = 2; cols2 ="deepskyblue" 
  # maxPeakDuration = "hrt";
  # interval_sec = 15
  
  
  # lagSec=6;
  # sgol_p = 2; sgol_n = 25;  weightMalus = 50;
  # match_threshold = 0.0; minSizeSec=4;
  # outputName = "PMBest_test";
  # correctionRangeSeconds = 0.5;
  # minPeakAmplitude = 0.05;
  # minPeakDuration = "hrt";
  # maxPeakDuration = "hrt"; #can be either a value in seconds, NA, or "hrt".
  # 
  
  
  debug = F
  # cols1 = 2; cols2 ="deepskyblue" 
  #   interval_sec = 15
    
    ###########################
    
    
    
    d1 = signal$s1
    d2 = signal$s2
    SR = frequency(signal)
    ransamp = lagSec * SR
    #setup weighting of similarity for distance.
    malus = weightMalus  / 100
    maxMalus = 1-weightMalus/100
    wx = (-ransamp):+ransamp
    weights_val = rangeRescale(dnorm(wx, mean=0, sd=1000),0,malus ) + maxMalus
    
    
    ### peaks-valleys detection ########
    pv1 = peakFinder(d1, sgol_p, sgol_n, mode = "b", correctionRangeSeconds, minPeakAmplitude)
    pv2 = peakFinder(d2, sgol_p, sgol_n, mode = "b", correctionRangeSeconds, minPeakAmplitude)
    
    
    # ai=355
    # xbest = signal$PMdev$xBest
    # ab = (pv1$samples[ai]-500):(pv1$samples[ai]+500)
    # rs1 = rangeRescale(signal$s1,0,1)
    # rs2 = rangeRescale(signal$s2,0,1)
    # plot(time(rs1)[ab],rs1[ab],t="l", ylim=c(0,1), main=session)
    # lines(time(rs2)[ab], rs2[ab],col=2)
    # points(time(rs1)[pv1$samples], rs1[pv1$samples],pch=2)
    # points(time(rs2)[pv2$samples], rs2[pv2$samples], col=2,pch=2)
    # 
    # points(time(rs1)[pv1$samples[pv1$class=="v"]], rs1[pv1$samples[pv1$class=="v"]])
    # points(time(rs2)[pv2$samples[pv2$class=="v"]], rs2[pv2$samples[pv2$class=="v"]], col=2)
    # segments(
    #   x0 = time(rs1)[xbest$s1],
    #   x1 = time(rs2)[xbest$s2],
    #   y0 = rs1[xbest$s1],
    #   y1 = rs2[xbest$s2]
    # )
    # rm(xbest)
    #####################################
    
    #' find triangles
    #' good triangles start with valleys. if not delete first element
    #' from bool and from all list elements of pv1 and pv2
    if(pv1$class[1] == "p"){
      pv1$bool[pv1$samples[1]] = FALSE
      for(k in 2:length(pv1)){
        pv1[[k]] = pv1[[k]][-1]
      }
    }
    if(pv2$class[1] == "p"){
      pv2$bool[pv2$samples[1]] = FALSE
      for(k in 2:length(pv2)){
        pv2[[k]] = pv2[[k]][-1]
      }
    }
    
    #series must end with valley. delete last peaks
    n1 = length(pv1$samples)
    if(pv1$class[n1] == "p"){
      pv1$bool[pv1$samples[n1]] = FALSE
      for(k in 2:length(pv1)){
        pv1[[k]] = pv1[[k]][-n1]
      }
    }
    n2 = length(pv2$samples)
    if(pv2$class[n2] == "p"){
      pv2$bool[pv2$samples[n2]] = FALSE
      for(k in 2:length(pv2)){
        pv2[[k]] = pv2[[k]][-n2]
      }
    }
    
    #just the peaks
    pp1 = data.frame(
      #the index must be rebuilt due to elimination of first last elements
      "index"   = (1:length(pv1$sample))[pv1$class=="p"],
      "samples" = pv1$samples[pv1$class=="p"],
      "time"    = pv1$time[pv1$class=="p"],
      "y"       = pv1$y[pv1$class=="p"],
      "amp"     = pv1$amp[pv1$class=="p"],
      "start_sample" = NA,
      "end_sample"   = NA
    )
    pp2 = data.frame(
      #the index must be rebuilt due to elimination of first last elements
      "index"   = (1:length(pv2$sample))[pv2$class=="p"],
      "samples" = pv2$samples[pv2$class=="p"],
      "time"    = pv2$time[pv2$class=="p"],
      "y"       = pv2$y[pv2$class=="p"],
      "amp"     = pv2$amp[pv2$class=="p"],
      "start_sample" = NA,
      "end_sample"   = NA
    )
    
    #find the maximum amplitudes for later scaling
    max1 = max(pp1$amp)
    max2 = max(pp2$amp)
    
    # M now ONLY contains p-peaks!
    M = matrix(NA ,nrow=length(pv1$samples),ncol=length(pv2$samples))
    
    # i=61
    
    #itera solo sui picchi p
    for(i in seq_along(pp1$sample)){
      
      # i e j sono gli iteratori riferiti alla tabella dei picchi
      # I e J sono gli iteratori riferiti a pv1 e pv2
      
      I = pp1$index[i]
      
      #trova il range attorno al picco I in cui cercare la lag
      search_range = (pp1$samples[i]-ransamp):(pp1$samples[i]+ransamp)
      
      # indice dei picchi in pp2 compresi in search_range
      matches_j = which(pp2$samples %in% search_range ) 
      
      # corrispondente indice rispetto a pv2
      matches_J = pp2$index[matches_j] 
      
      
      nMatch = length(matches_j)
      
      
      DEBUG = FALSE
      if(DEBUG ){
        oldPar = par()
        #create layout
        lrows = 3 #ceiling(nMatch/3)
        lcols = nMatch+1 #min(nMatch*2,2)
        mata = c(1,1,1, 1:(nMatch*3)+1)
        while(length(mata)%%3!=0){mata=c(mata,max(mata)+1)}
        layout(matrix(mata,byrow = T, ncol = 3 ))
        ab = max(1,pp1$samples[i]-ransamp-60):min(pp1$samples[i]+ransamp+60,length(d1))
        
        yfix1 = min(d1[ab]);yfix2 = min(d2[ab])
        myYlim = range(c(d1[ab],d2[ab]))
        myYlim = myYlim[2]-myYlim[1]
        par (mar=c(0,2,0,0))
        plot( time(d1)[ab], d1[ab]-yfix1,t="l",
              xlim = c(pv1$time[I]-lagSec-5 , pv1$time[I]+lagSec+5),
              ylim = c(0,myYlim), col=cols1)
        lines(time(d2)[ab], d2[ab]-yfix2,col=cols2)
        
        points(pv1$time,pv1$y-yfix1, col=cols1)
        points(pv1$time[I],pv1$y[I]-yfix1, pch=4, cex=3, col=cols1)
        
        points(pv2$time[pv2$samples %in% ab],pv2$y[pv2$samples %in% ab]-yfix2, col=cols2)
        points(pv2$time[matches_J],pv2$y[matches_J]-yfix2, pch=4, cex=3, col=cols2)
        
        
        arrows(pv1$time[I]-lagSec,pv1$y[I]-yfix1,pv1$time[I]+lagSec,code = 3)
        abline(v=c(pv1$time[I]-lagSec,pv1$time[I]+lagSec),lty=3)
        # suppressWarnings(par(oldPar));layout(matrix(1))
      }
      
      
      if(nMatch>0) {
        f=1
        for(f in 1:nMatch){ ###per ogni match
          
          #J è l'indice di pv2 a cui corrisponde il picco f
          J = matches_J[f]
          
          #j è l'indice di pp2 a cui corrisponde il picco f
          j = matches_j[f]
          
          lagv = pp2$samples[j] - pp1$samples[i] #distanza fra i due picchi, in samples.
          # Valori negativi indicano rosso anticipa blu, o s2 segue s1
          #         _/|_
          # s1_/\__/    \___ rosso
          #           _/|_
          # s2___/\__/    \_ blu
          
          ## MODE:3 
          ## applica la lag, poi intepola con un fattore medio tra v-p e p-v 
          ##
          ## es: v-p ratio: 2, p-v ratio: 4, mean ratio: 3
          ## d1:          v--p----v
          ## d2:        v----p----------------v
          ## d1':    v-------p-------------v
          ## d2':    v----p----------------v 
          ## keep:   |---------------------|
          ## stessa lunghezza (laggato) dell'altro picco.
          
          
          a1 = pv1$samples[I-1] 
          p1 = pv1$samples[I] 
          
          a2 = pv2$samples[J-1] 
          p2 = pv2$samples[J]
          
          # if decline is longer than half recovery time (hrt), use hrt
          # = amplitude of next valley is lower than amp/2
          if(maxPeakDuration == "hrt") {
            hrt_level = pv1$y[I] - pv1$amp[I]/2
            if(pv1$y[I+1]>hrt_level){
              b1 = pv1$samples[I+1]
            } else {
              #find sample corresponding to hrt_level
              decline = pv1$samples[I]:pv1$samples[I+1]
              b1 = which.min(abs(d1[decline]-hrt_level)) + pv1$samples[I]
            }
            if(DEBUG){
              abline(h=hrt_level-yfix1,lty=2,col=cols1)
              arrows(pv1$time[I-1],pv1$y[I-1]-yfix1,pv1$time[I-1],pv1$y[I] -yfix1,col=cols1)
              arrows(pv1$time[I],pv1$y[I]-yfix1,pv1$time[I], hrt_level-yfix1,col=cols1)
              points(time(d1)[b1],d1[b1]-yfix1,pch=13,cex=3,col=cols1)
            }
            
            #same for p2
            hrt_level = pv2$y[J] - pv2$amp[J]/2
            if(pv2$y[J+1]>hrt_level){
              b2 = pv2$samples[J+1]
            } else {
              #find sample corresponding to hrt_level
              decline = pv2$samples[J]:pv2$samples[J+1]
              b2 = which.min(abs(d2[decline]-hrt_level)) + pv2$samples[J]
            }
            if(DEBUG){
              abline(h=hrt_level-yfix2,lty=2,col=cols2)
              arrows(pv2$time[J-1],pv2$y[J-1]-yfix2,pv2$time[J-1],pv2$y[J] -yfix2,col=cols2)
              arrows(pv2$time[J],pv2$y[J]-yfix2,pv2$time[J], hrt_level-yfix2,col=cols2)
              points(time(d2)[b2],d2[b2]-yfix2,pch=13,cex=3,col=cols2)
            }
          } else if(!is.na(maxPeakDuration) && maxPeakDuration>0) {
            # maxPeakDuration è specificato e diverso da zero
            # controlla se c'è una valley prima, altrimenti tieni quel valore
            
            peakD1 = pv1$samples[I] + maxPeakDuration*SR
            b1 = min(pv1$samples[I+1], peakD1)
            
            peakD2 = pv2$samples[J] + maxPeakDuration*SR
            b2 = min(pv2$samples[J+1], peakD2)
            
            
          } else {
            #tieni tutto fino alla prossima valley
            b1 = pv1$samples[I+1]
            b2 = pv2$samples[J+1]
          }
          
          #' Interpola entrambi dalla v a sinistra fino a p e da p fino a hrt/v a destra
          l_ap = max((p1-a1),(p2-a2))+1
          l_pb = max((b1-p1),(b2-p2))+1
          
          ap1 = approx(x = 1:length(a1:p1), y = d1[a1:p1], xout=seq(1,length(a1:p1), length.out=l_ap))$y
          ap2 = approx(x = 1:length(a2:p2), y = d2[a2:p2], xout=seq(1,length(a2:p2), length.out=l_ap))$y
          pb1 = approx(x = 1:length(p1:b1), y = d1[p1:b1], xout=seq(1,length(p1:b1), length.out=l_pb))$y
          pb2 = approx(x = 1:length(p2:b2), y = d2[p2:b2], xout=seq(1,length(p2:b2), length.out=l_pb))$y
          
          toCor1 = c(ap1,pb1[-1])
          toCor2 = c(ap2,pb2[-1])
          
          #salva i valori di start/end utilizzati
          pp1$start_sample[i] = a1
          pp1$end_sample[i]   = b1
          pp2$start_sample[j] = a2
          pp2$end_sample[j]   = b2
          
          
          # stop("continue here")
          
          #adesso devi trovare una funzione per determinare la similarità fra i tuoi triangoli
          #idealmente dovrebbe già essere il valore di correlazione finale di AMICo
          #però devi riuscire a ritirare in ballo minSec cioè che micropicchi non si considerano da soli
          # e capire cosa fare per tutte le parti tra hrt e la valley successiva
          
          #abs(a-b) è la differenza tra le aree, non è male.
          #usare le derivate o no?
          #sto provando a normalizzare per i valori dei picchi individuali. sembra interessante!
          
          
          
          
          res1 = peakCor ( toCor1, toCor2, max1, max2,diff=F )
          res2 = crazyGold(toCor1, toCor2)
          res3 = cor(      toCor1, toCor2)
          res4 = normCor ( toCor1, toCor2, max1, max2)
          
          DEBUG= FALSE
          if(DEBUG){
            yfix1 = min(d1[a1:b1]) 
            yfix2 = min(d2[a2:b2])
            newx1 = (p1-l_ap+1):(p1+l_pb-1)
            newx2 = (p2-l_ap+1):(p2+l_pb-1)
            if(any(newx1<0)) newx1 = newx1 - min(newx1)+1
            if(any(newx1<0)) newx2 = newx2 - min(newx2)+1
            
            ymax =  max(d1[a1:b1]-yfix1,
                        d2[a2:b2]-yfix2)
            
            #original series (lagged)
            plot (time(d1)[a1:b1],      d1[a1:b1]-yfix1,t="l",col=cols1, ylim=c(0,ymax),xlim=range(time(d1)[newx1]))
            lines(time(d2)[a2:b2]-lagv/SR, d2[a2:b2]-yfix2,t="l",col=cols2)
            #peaks
            points(time(d1)[p1]     ,d1[p1]-yfix1,pch=17,cex=3,col=cols1)
            points(time(d2)[p2]-lagv/SR,d2[p2]-yfix2,pch=17,cex=3,col=cols2)
            
            # valori trasformati
            lines(time(d1)[newx1],      toCor1-yfix1,t="l",col=cols1,lwd=2,lty=2)
            lines(time(d2)[newx2]-lagv/SR, toCor2-yfix2,t="l",col=cols2,lwd=2,lty=2)
            
            #similarity:
            center = time(d1)[trunc(a1+ (b1-a1)/2)]
            text(center,
                 ymax/3,
                 paste0("peakCor:", round(res1,3), "\r\n",
                        "goldCor:", round(res2,3), "\r\n",
                        "Pearson:", round(res3,3), "\r\n",
                        "normCor:", round(res4,3), "\r\n"))
            
          }
          
          
          thisCor = res4
          ## WEIGHT CORRELATION
          
          # weight for distance
          weightCor = thisCor * weights_val[lagv+ransamp+1]
          
          # weight for stretching
          #' a logistic weight is applied according to the ratio
          #' of the longer and the shorter peaks
          #' if the longer is twice the shorter the malus is 0.7
          #' if it's 4 time, the malus is 0.5
          #' see them all: 
          # x = seq(1,10,by=0.1)
          # y =1/(1 + (x/350)^3.18)^(25*10^5)
          # plot(x,y );points(1:10,y[seq(1,100,by=10)],pch=16);y;abline(v=4,col=2)
          
          l2 = b2-a2
          l1 = b1-a1
          x = max(l2,l1)/min(l2,l1)
          y =1/(1 + (x/350)^3.18)^(25*10^5)
          weightCor = weightCor * y
          
          
          ## POPULATE FINAL MATRIX
          M[i,j] = weightCor
          
        }
      }
      
    }
    
    
    
    
    
    M[is.na(M)] = 0
    ## ignora le correlazioni negative. Non sono buoni match cmq! (?)
    M[M<0] = 0
    ## ignora le correlazioni sotto un certo threshold:
    M[M<match_threshold] = 0
    if(DEBUG){
      par(mar=c(2,2,0,0))
      layout(matrix(1))
      colors= c("#ffffff","#999999",  "#23c647")
      colfunc <- grDevices::colorRampPalette(colors=colors, bias=1)
      legendSteps =19
      colz = c("#ffffff",colfunc(legendSteps))
      iCol =round(rangeRescale(c(t(M)),1,20),0)
      nc = ncol(M)
      nr = nrow(M)
      
      plot(-1000, xlim=c(0,nc/8),ylim=c(0,nr/6),xaxs="i", yaxs="i")
      xs = 0:nc
      ys = 0:nr
      rect(
        xleft = xs[1:(length(xs)-1)],xright =xs[2:length(xs)],
        ybottom =  rep(ys[1:(length(ys)-1)], each = ncol(M)),
        ytop    =  rep(ys[2:length(ys)], each = ncol(M)),
        col=colz[iCol], border=NA
        
      )
      suppressWarnings(par(oldPar))
      
    }
    #### black magic starts here
    ## SORTING APPROACH n2: try all combinations
    ## https://cs.stackexchange.com/questions/91502
    
    #' REGOLE: ciascuna riga e ciascuna colonna possono essere matchate una sola volta
    #' inoltre non è possibile incrociare i match, se r2-c3, r3-c1 non è ammesso.
    #' Qua l'algoritmo crea una copia vuota di M, e comincia a popolarla iterativamente
    #' col valore più alto tra:
    #'   - la cella nella stessa riga ma colonna precedente
    #'   - la cella nella stessa colonna ma riga precedente
    #'   - la cella stessa + la cella diagonale precedente
    #'
    #' Questo permette di ottenere per ciascun match-riga colonna il valore di
    #' similarità cumulativo di quella scelta più la diagonale precedente, per trovare
    #' la soluzione che massimizza la similarità globale, e non del singolo picco.
    
    # new empty matrix
    A = matrix(rep(NA,length(M)), nrow = nrow(M) )
    best=list(row=rep(0,nrow(M)),col=rep(0,nrow(M)),similarity=rep(0,nrow(M)))
    #per tutte le celle di M, per riga:
    
    for (i in 1:nrow(M)){
      for(j in 1:ncol(M)){
        A[i,j] = max( max(A[i-1,j],0), #same row previous cell
                      max(A[i,j-1],0), #previous row, same cell
                      max(A[i-1,j-1],0)+M[i,j],#previous diagonal + current cell
                      na.rm = T)
      }
    }
    
    if(DEBUG){
      
      # par(mar=c(2,2,0,0))
      layout(matrix(1))
      colors= c("#ffffff","#999999",  "#23c647")
      colfunc <- grDevices::colorRampPalette(colors=colors, bias=1)
      legendSteps =19
      colz = c("#ffffff",colfunc(legendSteps))
      iCol =round(rangeRescale(c(t(A)),1,20),0)
      nc = ncol(A)
      nr = nrow(A)
      
      plot(-1000, xlim=c(0,nc),ylim=c(0,nr),xaxs="i", yaxs="i")
      xs = 0:nc
      ys = 0:nr
      rect(
        xleft = xs[1:(length(xs)-1)],xright =xs[2:length(xs)],
        ybottom =  rep(ys[1:(length(ys)-1)], each = ncol(A)),
        ytop    =  rep(ys[2:length(ys)], each = ncol(A)),
        col=colz[iCol], border=NA
        
      )
      suppressWarnings(par(oldPar))
      
    }
    
    #' Ora l'algoritmo procede al contrario, dall'angolo alto dx, ovvero dai valori
    #' più alti della matrice. Per ciascuna riga e colonna trova il massimo
    #' (che implica che la somma dei valori ammessi precedenti sia massimo)
    #' salva quel valore, ed elimina tutto su quella riga e quella colonna (perché
    #' sono già stati scartati)
    
    A2 = A
    nr = nrow(A2)
    nc = ncol(A2)
    counter = 0
    if(DEBUG){
      M2 = M
      scalemax = max(A2)
    }
    while(sum(A2)){
      counter = counter +1
      x = which(A2 == max(A2), arr.ind = TRUE)[1,] #there might be more than 1 result. Choose the first.
      A2[x[1]:nr,] = 0
      A2[,x[2]:nc] = 0
      #' x[1] per pp1 e x[2] per pp2 rappresentano l'indice del picco selezionato
      best$row[counter]= x[1] #pp1 index
      best$col[counter]= x[2] #pp2 index
      best$similarity[counter] = M[x[1],x[2]] # similarity
      
      
      ## This code creates the frames of an animation ####################
      # if(DEBUG){
      # 
      # 
      #   png(paste0("dev_spam/AMICo matrix/",counter,".png"),width = 800, height = 800, units = "px", type="cairo" )
      #   
      #   par(mar=c(2,2,0,0))
      #   layout(matrix(1))
      #   colors= c("#ffffff","#999999",  "#23c647")
      #   colfunc <- grDevices::colorRampPalette(colors=colors, bias=1)
      #   legendSteps =19
      #   colz = c(rgb(1,1,1,0),colfunc(legendSteps))
      #   iCol =round(rangeRescale(c(t(M2)),1,20),0)
      #   nc = ncol(M2)
      #   nr = nrow(M2)
      #   
      #   plot(-1000, xlim=c(x[2]-20,x[2]+20),ylim=c(x[1]-20,x[1]+20),xaxs="i", yaxs="i")
      #   rect(x[2]-1, -10000, x[2]+20,x[1]+20,
      #        col="#333333",border=NA)
      #   rect(x[2]-20, x[1]-1, x[2]+20,x[1]+20,
      #        col="#333333",border=NA)
      # 
      #   
      #   xs = 0:nc
      #   ys = 0:nr
      #   rect(
      #     xleft = xs[1:(length(xs)-1)],xright =xs[2:length(xs)],
      #     ybottom =  rep(ys[1:(length(ys)-1)], each = ncol(M2)),
      #     ytop    =  rep(ys[2:length(ys)], each = ncol(M2)),
      #     col=colz[iCol], border=NA
      #     
      #   )
      #   rect(x[2]-1, -10000, x[2]+20,x[1]+20,
      #        col="#333333",border=NA,angle=0, density = 20)
      #   rect(x[2]-20, x[1]-1, x[2]+20,x[1]+20,
      #        col="#333333",border=NA,angle=0, density = 20)
      #   rect(
      #     x[2]-1, x[1]-1, x[2],x[1],
      #     col=NA, border=2,lwd=4
      #   )
      #   dev.off()
      # 
      #   
      # } ####################
    }      
    
    
    
    xbest = data.frame(best)
    xbest = xbest[order(xbest$row),]
    xbest = xbest[xbest$col!=0 & xbest$row !=0,]
    xbest$lag =   pp2$samples[xbest$col]-
      pp1$samples[xbest$row]
    xbest$a1 =    pp1$start_sample[xbest$row]
    xbest$p1 =    pp1$samples[xbest$row]
    xbest$b1 =    pp1$end_sample[xbest$row]
    
    xbest$a2 =    pp2$start_sample[xbest$col]
    xbest$p2 =    pp2$samples[xbest$col]
    xbest$b2 =    pp2$end_sample[xbest$col]
    xbestrealn = nrow(xbest)
    
    if(xbestrealn < 2) {
      mess = paste("In session", UID(signal), "no matches were found. Please check the raw data and filtering pipeline.\n")
      # cat(mess, file=stdout())
      warning(mess)
      return(newAMICo(sync = NA,NA,xbest,args))
    } else {
      ############################################################################################
      #' now now now
      #' we know the similarity in proximity of the peaks. What about in between?
      #' Proposal [a]: we use the same logic! interpolate the inbetweens, normalize them
      #' and calculate the similarity. 
      
      # something special for the parts before the first peak
      # warning("to do")
      
      ampz1 = ampz2 = c()
      intsamp = interval_sec*SR
      minsize = minSizeSec*SR
      
      for(i in 2:nrow(xbest)){
        if(i==2 && nrow(xbest) != xbestrealn) stop("resetta xbest")
        #inbetween of series 1
        ba1 = xbest$b1[i-1]:xbest$a1[i]
        ba2 = xbest$b2[i-1]:xbest$a2[i]
        if (all(c(length(ba1),length(ba2)) == 1)){
          # non c'è niente da interpolare!
          # Nothing to look here, move on.
        } else if (all(c(length(ba1),length(ba2)) > minsize)){
          #ok qua c'è da fare! entrambe le due serie hanno del materiale in mezzo
          # stop()
          
          ampz1[i] = max(d1[ba1]-min(d1[ba1]))
          ampz2[i] = max(d2[ba2]-min(d2[ba2]))
          
          
          if(DEBUG){
            par(mfrow=c(1,1))
            
            plot(time(d1)[ba1],d1[ba1]-min(d1[ba1]), col=cols1,t="l",
                 ylim=c(0,max(c(d1[ba1],d2[ba2]))),
                 xlim=c(min(c(time(d1)[ba1],time(d2)[ba2])),max(c(time(d1)[ba1],time(d2)[ba2]))))
            lines(time(d2)[ba2],d2[ba2]-min(d2[ba2]), col=cols2)
          }
          #1 dividi per intervalli di max interval_sec creando n split (contando sul più breve)
          # if(any(time(d1)[ba1]>40*60+55)) stop()
          
          #duration of the shortest segment
          minl = min(length(ba1),length(ba2))
          #in how many splits can we split the shortest?
          splits = max(1,floor(minl/intsamp))
          
          #durations of the splits
          sdur1 = floor(length(ba1)/splits) #durata split s1
          sdur2 = floor(length(ba2)/splits) #durata split s2
          sdurmax = max(sdur1, sdur2)       #durata split finale
          s=1
          for(s in 1:splits){
            
            # i sample dello split corrente
            sba1 = ba1[(1:sdur1)+(s-1)*sdur1]
            sba2 = ba2[(1:sdur2)+(s-1)*sdur2]
            
            #2 interpola ciascuna finestra alla stessa lunghezza
            ib1 = approx(x = 1:sdur1, y = d1[sba1], xout=seq(1, sdur1, length.out=sdurmax))$y
            ib2 = approx(x = 1:sdur2, y = d2[sba2], xout=seq(1, sdur2, length.out=sdurmax))$y
            
            if(DEBUG){
              par(mfrow=c(2,1))
              
              plot(rangeRescale(d1[ba1]-min(d1[ba1]),0,1,0,max1),
                   main=round(cor(ib1, ib2),3), ylim=c(0,1))
              lines(1:sdurmax+(s-1)*sdurmax, rangeRescale(ib1-min(d1[ba1]),0,1,0,max1), col=cols1, lty=2,lwd=3)
              lines((1:sdur1)+(s-1)*sdur1, rangeRescale(d1[sba1]-min(d1[ba1]),0,1,0,max1))
              
              
              points(rangeRescale(d2[ba2]-min(d2[ba2]),0,1,0,max2))
              lines(1:sdurmax+(s-1)*sdurmax, rangeRescale(ib2-min(d2[ba2]),0,1,0,max2), col=cols2, lty=2,lwd=3)
              lines((1:sdur2)+(s-1)*sdur2, rangeRescale(d2[sba2]-min(d2[ba2]),0,1,0,max2))
              
              mtext(paste("norm cor:", round(normCor(ib1, ib2, max1, max2),3),
                          "simple:",round(sigCor(ib1, ib2,max1,max2),3)#,
                          # "mutin:",round(mutin(ib1, ib2,max1,max2),3)
              ))
              
            }
            
            #4 calcola similarity (nomrlizza + RMSD delle derivate)
            #dal momento che non ci sono picchi qua, il massimo teorico è 
            # minpeak delta
            # res = peakCor(ib1, ib2, max(ib1), max(ib2),diff=T)
            # res = crazyGold(ib1,ib2)
            res = normCor(ib1,ib2,max1,max2 )
            
            # print(crazyGold(ib1, ib2))
            
            # plot(ib2-min(ib2),ib1-min(ib1))
            # sum(residuals(lm(I(ib2-min(ib2))~I(ib1-min(ib1))))^2)
            
            # weight for distance
            weightCor = res * weights_val[abs(length(sba1)-length(sba2))+ransamp+1]
            
            # weight for stretching
            
            x = max(sdur1,sdur2)/min(sdur1,sdur2)
            y =1/(1 + (x/200)^3.18)^(25*10^5)
            weightCor = weightCor * y
            if(DEBUG){
              mtext(paste("weighted:",round(weightCor,3)) )
            }
            
            # #entrambe le due serie hanno del materiale in mezzo
            # #però almeno una è troppo corta. Mettiamo tutto a NA
            # if(min(length(sba1),length(sba2))<minsize*0.7){
            #   weightCor = NA
            # }
            
            
            #5 crea delle righe aggiuntive in xbest con questi valori
            p1 = sba1[round(sdur1/2)]
            p2 = sba2[round(sdur2/2)]
            xbest=rbind(xbest,data.frame(
              "row" = NA,
              "col" = NA,
              "similarity" = weightCor,
              "lag" = p2 - p1,
              "a1"  = sba1[1],
              "p1"  = p1,
              "b1"  = sba1[length(sba1)],
              "a2"  = sba2[1],
              "p2"  = p2,
              "b2"  = sba2[length(sba2)]
            ))
            
          }
          
        } else {
          #'  qua entrambe le serie sono più lunghe di 1, ma sono troppo corte per calcolare
          #' la sincro. Oppure una valley è usata 2 volte, mentre l'altro segnale ha della roba in mezzo
          #' 
          #'     /\  /\
          #' s1 /  \/  \
          #' 
          #'    /\     /\
          #' s2/  v\__/  \
          #' 
          #' In alternativa possiamo semplicemente stretchare i due punti del segnale
          #' con la roba in mezzo al loro midpoint:
          midpoint = round(mean(ba2))
          xbest$b2[i-1] = midpoint
          xbest$a2[i] = midpoint
          
          midpoint = round(mean(ba1))
          xbest$b1[i-1] = midpoint
          xbest$a1[i] = midpoint
        }
      } 
      
      # something special for the parts after the last peak
      # warning("to do")
      
      
      ############################################################################################
      
      #' Now xbest is done. Create the 2 interpolated series, sync & lag
      syncvec = lagvec = rep(NA,length(d1))
      xt = time(d1)
      xbest = xbest[order(xbest$a1 ),]
      for(i in 1:nrow(xbest)){
        #midpoint between s1 and s2
        # a = round(mean(xt[xbest[i,"a1"]], xt[xbest[i,"a2"]]))
        # b = round(mean(xt[xbest[i,"b1"]], xt[xbest[i,"b2"]]))
        xbest$a[i] = round(mean(c(xbest[i,"a1"], xbest[i,"a2"])))
        xbest$b[i] = round(mean(c(xbest[i,"b1"], xbest[i,"b2"])))
        xbest$ta[i]= xt[xbest$a[i]]
        xbest$tb[i]= xt[xbest$b[i]]
        a = round(xbest$ta[i])
        b = round(xbest$tb[i])
        syncvec[a:b] = xbest$similarity[i]
        lagvec[a:b]  = xbest$similarity[i]
      }
      

      
      #' XBEST guide:
      #' row: the peak number of s1 wich was matched
      #' col: the peak number of s2
      #' similarity: some similarity computation
      #' lag: the delta between the sample of p1 and p2
      #' a1, p1, b1: sample position of onset, peak, and end of a feature
      #' a2, p2, b2: the same for s2
      #' a, b: the average start and end between s1 and s2
      #' ta, tb: the time corresponding to a and b

      # finallyt instantiate new sync class object
      sync = rats(syncvec, start = start(d1),
                  frequency=frequency(signal), timeUnit="second",
                  valueUnit="0-1")
      lags = rats(lagvec,  start = start(d1),
                  frequency=frequency(signal), timeUnit="second",
                  valueUnit="seconds")
      # applica gli artefatti
      for(i in seq_len(nrow(signal$artefacts)) ){ 
        window(sync, signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
        window(lags, signal$artefacts[i,"start"], signal$artefacts[i,"end"] ) <- NA
      }
      
      return(newAMICo(sync,lags, xbest, args))
      
    }
    
}


peakCor = function(x,y,xmax,ymax, diff){
  #' this similarity index rescales the v-p-v triangle according to the 
  #' highest peak observed by that participant during the session
  #' then calculates the root mean square of the absolute differences
  DEBUG = FALSE
  ###DEBUG
  # stop("DEBUG")
  # x = ib1
  # y = ib2
  # xmax = minPeakDelta
  # ymax = minPeakDelta
  # diff=T
  ########
  if(length(x)!=length(y))stop("x and y must have same length")
  x = rangeRescale(x-min(x),0,1, 0, xmax)
  y = rangeRescale(y-min(y),0,1, 0, ymax)
  if(diff){
    x = diff(x)
    y = diff(y)
  }
  #'ora x e y sono normalizzati su una scala da 0 a 1 dove 1 è il picco
  #'più alto di quel soggetto.
  #'Dimostrazione: max(x)*xmax == pv1$amp[i]
  
  
  
  #' YET ACTUALLY
  #' the maximum different situation is that of a parabola with apex (M,1)
  #' compared to a flat zero line or parabola with apex (M,0)
  #' so the maximum possible difference is the first parabola area
  #' Area = 2/3 * 2M, where 2M equals the length of x
  #'
  
  #parabola
  l = length(x)
  y2 = -1/(l/2)^2 * (1:l - l/2)^2 +1
  if(diff){
    y2 = diff(y2)
  }
  max_p = sqrt(sum(y2^2)/length(y2))
  #' In this new crazy formula, I calculate the maximum using the parabola
  #' and normalize for that maximum to get back to a 0-1 value
  (sqrt(sum(y2^2)/l) - sqrt(sum(abs(x-y)^2)/l)) /  sqrt(sum(y2^2)/l)
  
  
  
  #' this still works Best:
  #' 1 - the Root-mean-square deviation of the 2 peaks
  #' max_p = 1
  
  res = (max_p-sqrt(sum((x-y)^2)/length(x))) / max_p
  
  if(DEBUG){
    plot (x, col=cols1, t="l",ylim=c(min(c(x,y,y2)),max(c(x,y,y2))), main=res)
    lines(y, col=cols2)
    segments(1:length(x), x,
             1:length(y), y)
    if(!is.null(y2)){
      lines(y2, lwd=2, col=3)
    }
  }
  res
  
}

crazyGold = function(x,y){
  DEBUG= FALSE
  if(length(x)!=length(y))stop("x and y must have same length")
  
  d1 = as.numeric(diff(x)); d2 = as.numeric(diff(y))
  d1 = rangeRescale(d1,0,1); d2 = rangeRescale(d2,0,1)
  
  if(DEBUG){
    # y = -x
    plot (d1, col=cols1, t="l",ylim=c(0,1))
    lines(d2, col=cols2)
    segments(1:length(d1), d1,
             1:length(d2), d2)
    # plot(1-abs(d1 - d2)^2)
    
  } 
  
  (1-mean(abs(d1 - d2)))^2
  
}


crazyCor = function(x,y){
  DEBUG= FALSE
  x = as.numeric(diff(x)); y = as.numeric(diff(y))
  
  if(length(x)!=length(y))stop("x and y must have same length")
  if(DEBUG){
    # y = -x
    plot (rangeRescale(x,0,1),t="l")
    lines(rangeRescale(y,0,1),col=2)
    segments(1:length(x),rangeRescale(x,0,1),
             1:length(y),rangeRescale(y,0,1))
  }
  m2d = 1-mean(abs(rangeRescale(x,0,1) - rangeRescale(y,0,1)))
  mean(m2d)^2
  
  1- sqrt(sum(abs(rangeRescale(x,0,1) - rangeRescale(y,0,1))^2))
  
}


normCor = function(x,y,xmax=max(x)-min(x),ymax=max(y)-min(x)){
  #https://www.youtube.com/watch?v=ngEC3sXeUb4
  #this is actually the cosine similarity
  if(length(x)!=length(y))stop("x and y must have same length")
  x = rangeRescale(x-min(x),0,1, 0, xmax)
  y = rangeRescale(y-min(y),0,1, 0, ymax)
  
  sum(x*y) / sqrt(sum(x^2)*sum(y^2))
  
}

# x = 1:10
# y = 10:1
# plot(x,t="l")
# lines(y, col=2)
# xmax = ymax = 10
# if(length(x)!=length(y))stop("x and y must have same length")
# x = rangeRescale(x-min(x),0,1, 0, xmax)
# y = rangeRescale(y-min(y),0,1, 0, ymax)
# 
# sum(x*y) / sqrt(sum(x^2)*sum(y^2))
# y = y-min(y)
# plot(x,t="l")
# lines(y, col=2)

sigCor = function(x,y,xmax=max(x)-min(x),ymax=max(y)-min(x)){
  #https://www.youtube.com/watch?v=ngEC3sXeUb4
  x = rangeRescale(x-min(x),0,1, 0, xmax)
  y = rangeRescale(y-min(y),0,1, 0, ymax)
  sum(x*y)/length(x)
  
}

mutin = function(x,y,xmax,ymax){
  if(!missing(xmax) && !missing(ymax)){
    x = rangeRescale(x-min(x),0,1, 0, xmax)
    y = rangeRescale(y-min(y),0,1, 0, ymax)
  }
  copent::copent(data.frame(x,y))/sqrt(abs(copent::entknn(x)*
                                             copent::entknn(y)))
}



