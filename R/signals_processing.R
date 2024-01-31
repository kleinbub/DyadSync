#' Finds peaks and throughs features in a time series
#'
#' @param x a rats, ts, or a numeric vector.
#' @param sgol_p smoothing filter order. If NA, the filtering is disabled.
#' @param sgol_n smoothing filter length (must be odd).
#' @param mode should the function return only 'peaks', 'valleys', or 'both'?
#' @param correctionRangeSeconds the half range in which the real maximum/minimum
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
  which.minmax = function(x, pv){if(pv=="p") which.max(x) else if(pv=="v") which.min(x)}
  
  if(missing(sgol_p)) stop("in peakfinder missing sgol_p")
  if(missing(sgol_n)) stop("in peakfinder missing sgol_n")
  if(missing(correctionRangeSeconds)) stop("in peakfinder missing correctionRangeSeconds")
  if(missing(minPeakAmplitude)) stop("in peakfinder missing minPeakAmplitude")
  
  SR = frequency(x)
  timeX = time(x)
  attributes(x) = NULL
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
  piks = timeX[piksam]
  
  if(length(piksam) == 0) warning("No peaks were found with the current settings!")
  
  list("bool" = pikboo,
       "samples" = piksam,
       "time" = piks,
       "class" = pv,
       "y" = x[piksam],
       "amp" = c(NA,diff(x[piksam])),
       "index" = is
  )
  
}