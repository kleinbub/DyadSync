##     ___                 _   ___ _
##    /   \_   _  __ _  __| | / __\ | __ _ ___ ___
##   / /\ / | | |/ _` |/ _` |/ /  | |/ _` / __/ __|
##  / /_//| |_| | (_| | (_| / /___| | (_| \__ \__ \
## /___,'  \__, |\__,_|\__,_\____/|_|\__,_|___/___/
##         |___/
############################################################################################################  
## DyadRandom.R
# bootstrap and permutation tools
#
############################################################################################################ 
## ToDo
# -
#
############################################################################################################ 

#getCombo calcola le possibili combinazioni di segnali, escludendo solo
#i segnali realmente appaiati. Non tiene in considerazione il ruolo,
#quindi nuove diadi composte dal segnale di due pazienti sono possibili.

getCombo = function(n_dyads, max_combo=10000){
  elements = c(paste0("s1",":",1:n_dyads),paste0("s2",":",1:n_dyads))
  combo = combn(elements,2)
  comborem = sapply(1:ncol(combo), function(i){strsplit(combo[1,i],":")[[1]][2] != strsplit(combo[2,i],":")[[1]][2] })
  combo = combo[,comborem]
  #if(verbose)print(combo)
  if (ncol(combo) <= max_combo) { #keep all the combinations
    com = data.frame(t(combo),stringsAsFactors=F)
    #com =com[sample(nrow(com)),] #shuffle
  } else { 
    zamp = sample(1:ncol(combo), size=max_combo ) #keep only max_combo combinations
    com = data.frame(t(combo[,zamp]),stringsAsFactors=F)
  }
  rownames(com)=1:nrow(com)
  cat(paste("\r\n",nrow(com),"drawn out of",nCombo(n_dyads),"possible combinations"))
  print(table(as.vector(as.matrix(com))))
  par(las=2,oma=c(2,0,0,0))
  plot(table(as.vector(as.matrix(com))),main="Signal instances in sample", ylab = "number of instances")
  com$ranx_session = sapply(com$X1, function(x) as.numeric(strsplit(x,":")[[1]][2]))
  com$ranx_role = sapply(com$X1, function(x) strsplit(x,":")[[1]][1])
  com$rany_session = sapply(com$X2, function(x) as.numeric(strsplit(x,":")[[1]][2]))
  com$rany_role = sapply(com$X2, function(x) strsplit(x,":")[[1]][1])
  com = com[ , 3:6]
  ##print(table(c(as.numeric(com$ranx_session),as.numeric(com$rany_session))))
  return(com)
}

#' Calculate the number of possible pseudo-dyads
#'
#' @param n_dyads 
#' @param max_combo 
#'
#' @return
#' @export
#'
#' @examples
nCombo = function(n_dyads, max_combo=10000){
  n = n_dyads*2
  r=2
  cmb = factorial(n)/(factorial(r)*factorial(n - r)) -n_dyads
  if(max_combo<cmb) cmb = max_combo
  return(cmb)
}


#############################################
## Bootstrap and permutation tools
##
## the following functions take real dyads data and outputs randomized dyads
## These are useful to compare observed synchrony with random synchrony.
## dyadComb uses a combination approach (bootstrap)
## dyadPerm uses permutations, which allow for exact testing.
## dyadShift instead of shuffling between dyads, shifts and reverse the data inside each dyad.

#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
dyadComb = function(x,...) {UseMethod("dyadComb")}

#' @export
dyadComb.DyadExperiment = function(exper, n, signals="all", verbose=F){
  
  #Questa funzione accetta in ingresso un esperimento, e genera un nuovo esperimento con n "sedute"
  #random, in cui vengono scombinati gli appaiamenti paziente/terapeuta
  
  
  #debug
  # exper = lr
  # n =10
  # signals="all"
  # verbose = T
  
  if(sum(signals == "all")){
    signals = names(exper[[1]])
    nSignals = length(signals)
  } else {
    if(is.vector(signals)) {
      nSignals = length(signals)
      exper = selectSignals(exper,signals)
    } else stop("signals is not well specified")
  }
  if(sum(vapply(exper[[1]],is.DyadSignal,c(T))) != nSignals) stop ("Not all specified signals are DyadSignal objects")
  
  
  ncom = nCombo(length(exper),n)
  
  #create a different reshuffling for each signal
  comboList = list()
  for(i in 1:length(signals)){
    comboList[[i]]=getCombo(length(exper),n ) 
  }
  names(comboList)=signals
  if(verbose) print(comboList)
  
  #do the magic
  sessionList = list()
  for(i in 1:ncom){
    #cat("\r\n genero la seduta posticcia ",i)
    signalList = list()
    for(j in 1:nSignals){
      if(is.DyadSignal(exper[[1]][[signals[j]]])){
        #estraggo lo stream raw per la sessione ed il ruolo identificati in combolist
        jSampRate = frequency(exper[[1]][[signals[j]]])
        olds1name = s1Name(exper[[ comboList[[j]][i,"ranx_session"] ]])
        olds2name = s2Name(exper[[ comboList[[j]][i,"ranx_session"] ]])
        
        s1raw = as.numeric(exper[[ comboList[[j]][i,"ranx_session"] ]][[signals[j]]][[comboList[[j]][i,"ranx_role"]]])
        s2raw = as.numeric(exper[[ comboList[[j]][i,"rany_session"] ]][[signals[j]]][[comboList[[j]][i,"rany_role"]]])
        # tieni la lunghezza del più corto
        s1raw = rats(s1raw[1:min(length(s1raw),length(s2raw))], frequency = jSampRate, timeUnit="second")
        s2raw = rats(s2raw[1:min(length(s1raw),length(s2raw))], frequency = jSampRate, timeUnit="second")
        #ricostruisci il dyadsignal

        newSignal = DyadSignal(name= signals[j],s1 = s1raw, s2 = s2raw, SR = jSampRate,
                               s1Name= paste0(comboList[[j]][i,"ranx_session"], if(comboList[[j]][i,"ranx_role"] =="s1") olds1name else olds2name ),
                               s2Name =paste0(comboList[[j]][i,"rany_session"], if(comboList[[j]][i,"rany_role"] =="s1") olds1name else olds2name ),
                               sessionId = i, dyadId = "random", groupId="random")
        signalList[[j]] = newSignal
        
      } else if(is.DyadCategory(exper[[1]][[signals[j]]])) {
        ##code for dyadcategories
        
      } else {
        warning("Unrecognized object was skipped in session",i, "at index", j)
        
      }
      
    }
    names(signalList) = signals
    
    sessionList[[i]] = DyadSession(
      groupId = "randomGroup",
      sessionId = paste0("randomSession_",i),
      dyadId = "randomDyad",
      signalList = signalList, s1Name ="random", s2Name = "random",fileName = "random")
  }
  names(sessionList) = paste0("shuffle_",1:ncom)
  nexp = DyadExperiment(paste0(n,"_combinations"),sessionList)
  #cat('\r\n',nrow(com),"/",length(combo)/2,"possible combinations where randomly selected")
  return(nexp)
}

########################CCF COMPARE
##This function plots the original and the random experiments
ccfCompare = function(EXP, random, signal="SC", lag="bestCCF", FUN ="marciIndex"){
  FUNname = FUN
  FUN <- match.fun(FUN)
  #print("Result is a matrix with each row rapresenting a random dyad, each column a window")
  resList = lapply(EXP, function(ses){ses[[signal]]$ccf$ccfmat[[lag]]})
  res = do.call("rbind",resList)
  colnames(res) = paste0("w",1:ncol(res))
  
  if(!is.null(random)){
    ranList = lapply(random, function(ses){ses[[signal]]$ccf$ccfmat[[lag]]})
    ras = do.call("rbind",ranList)
    colnames(ras) = paste0("w",1:ncol(ras))
  }
  
  #applica la funzione su ciascuna seduta
  marRes = apply(res,1,FUN)
  marRas = apply(ras,1,FUN)
  
  # cohend(marRes,marRas)
  # cohend(as.vector(res),as.vector(ras))
  
  if(FUNname == "marciIndex") FUNname = "Marci index"
  else FUNname = paste(FUNname, "x-cor")
  
  par(mfrow=c(1,2),oma=c(0,0,2,0))
  boxplot(list("Real"=marRes,"Random"=marRas), lty=c(1,1),ylab= FUNname,
          main="")
  
  #abline(h=c(0,0.5,-0.5))
  
  plot(density((marRes))$x,density((marRes))$y,ylim=c(0,max(density(marRes)$y,density(marRas)$y )),
       main= "", xlab=paste("N (real) =", length(EXP),"| N (random) =", length(random)), type="l",ylab= paste(FUNname,"density"))
  
  
  lines(density((marRas)),lty=3)
  legend("topright",legend=c("Real","Random"), lty=c(1,3))
  title(main = "Real and random dyads", outer= T,line=-1)
  title(main =paste("Signal:",signal,"| lag:",lag), 
        outer= T, cex.main=0.9,line = -2)
  title(main= paste("Cohen's d:",round(cohend(marRes,marRas),2),
                    "| t:",round(t.test(marRes,marRas)$statistic,2), 
                    " p-val:",pvalFormat(t.test(marRes,marRas)$p.val)),
        cex.main=0.8,line=-3,outer=T)
  
  #print(quantile(res))
  #return(res)
}
