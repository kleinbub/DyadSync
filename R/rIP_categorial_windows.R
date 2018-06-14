##stato del file:

#il file contiene 4 funzioni:
#----------------------------------
#catSummary
# questa è piuttosto ok
# -è inefficiente l'estrazione dello stream ad ogni livello categoriale -_-"
# -per carità dividi il motore del "merged" (o output full) da quello delle statistiche
# -in catSummary.dyadExperiment metti un check sulle categorie con differenza < 3 caratteri (su stringhe >3 caratteri)
#
#----------------------------------
#catRandom
# per funzionare funziona. però aggregare mille mila campionamenti 
# alla fine è abbastanza uguale ad usare la media/mediana del campione originale
# ToDo: Sta roba sarebbe più utile per calcolare intervalli di confidenza sulle misure

#----------------------------------
#plotCatLag
# plotta nello stesso boxplot confronti di lag e ccf tra categorie
# funziona bene ma plotta una session alla volta
# ToDo:
# -permetti il sort delle categorie
# -plot longitudinale e/o confronto diretto tra due casi
#
#----------------------------------
#plotCatStream
# siccome plotCatLag faceva sanguinare gli occhi e contorcere le menti di tutti
# plotCatStream è una funzione generica che mostra UN solo stream (es: lag o ccf)
# e non più le due informazioni insieme. Inoltre è più agile.
# ToDo:
#   -sviluppare plotCatStream.DyadSession
#   -sviluppare il grafico di comparison fra più casi in plotCatStream.DyadExperiment
#
#----------------------------------


#marci index ottimizzato. (taglia i valori infiniti a +- 10)
#per far na roba fatta ben si dovrebbe vedere la distribuzione casuale di marciIndex e vedere come e dove tagliare
marciIndex = function(a){
  if(is.list(a)){ ##se 'a' è una lista, applica l'indice a ogni elemento della lista
    matrix(unlist(lapply(a, function(iFile){apply(iFile,2,marciIndex)})),nrow=length(a),byrow=T)
  } else {
    x=log(sum(a[a >0]) / abs(sum(a[a<0]))) #formula dell'indice
    if(!is.na(x)){
      if(x>=10)x=10 else if(x <= -10)x=-10
    }
    
    x
  }
}

#catsummary ora dovrebbe essere più dignitosa, calcolando
#media e dispersione di tutte le finestre di ciascuna categoria INSIEME
# seprando solo per sessione (e non più per occorrenza)

catSummary <- function(x, signal="RRmean", streamKey="bestCCF", category="PACS", column="CATEGORIA",return.type=c("summary","full")) {
  if(!is.DyadExperiment(x) & !is.DyadSession(x))stop("only experiments or sessions get butter")  
  UseMethod("catSummary",x)
}

catSummary.DyadExperiment=function(experiment, signal="SC", streamKey="bestCCF", category="PACS", column="CATEGORIA",return.type=c("summary","full")) {
  warning("This is bugged and not yet supported. Use catSummary only for individual sessions")
  res = list()
  for(session in experiment){
    res[[paste(attr(session,"sessionId"),attr(session,"dyadId"),sep="_")]] = c(res,catSummary(session,signal,streamKey,category,column))
  }
  
    return(res)
}

catSummary.DyadSession=function(session, signal="SC", streamKey="bestCCF", category="PACS", column="CATEGORIA",return.type=c("summary","full")) {
  #####debug
  # session = full[[1]]
  # signal="SC"
  # streamKey="bestCCF"
  # category="PACS"
  # column="CATEGORIA"
  # return.type="summary"
  
  ###
  ### CatSummary per ciascun livello categoriale estrae lo stream streamKey, lo taglia con i minutaggi di ciascuna istanza di quel livello
  ### poi incolla insieme tutti i tagli e restituisce un vettore con tutti i valori incollati
 
  return.type=match.arg(return.type)
  #§DEBUG:
  #streamKey = stream = "bestLag"
  #dots= list(...)
  #print(dots)
  ##step 1:individuare TUTTI i livelli possibili in tutto l'esperimento
  # allCat = factor(as.character(unlist(sapply(full, function(session){as.character(session$categ[[category]][,column])}))))
  allLevels = levels(session$categ[[category]][[column]])
  # full2= lapply(full,function(session){
  #   session$categ[[category]][,column] = factor(as.character(session$categ[[category]][,column]),levels=allLevels)
  #   session
  # })
  # attributes(full2) =attributes(full)
  
  res = list()
  #per ciascun livello categoriale, ritaglia le finestre e calcola gli indici
  for(a in allLevels){
    sessions= data.frame("session"=0,"dyadId"=0,"mean"=0,"sd"=0,"median"=0,"marci"=0,"duration"=0,
                         "min_duration"=0,"max_duration"=0, "count"=0,"quart1"=0,"quart3"=0,
                         "min"=0,"max"=0)
    i=1
    #for(session in full2){
    #session = full2[[i]]
    print("-------------------------------------------------------")
    print(a)
    #print(paste(attr(session,"sessionId"),attr(session,"dyadId")))
    
    #this attempt to identify the stream by strings instead of by precise column name... WHY?!
    capture.output({
      if(grepl("pat",streamKey,ignore.case = T) || grepl("paz",streamKey,ignore.case = T)){
        stream = session$signals[[signal]]$patient
        #print("grepl PATIENT")
      } else if (grepl("ter",streamKey,ignore.case = T)|| grepl("ther",streamKey,ignore.case = T) || grepl("cli",streamKey,ignore.case = T)){
        stream = session$signals[[signal]]$clinician
        #print("grepl CLINICIAN")
      } else {stream = getCCF(session$signals[[signal]], streamKey)
      #print(paste("grepl",streamKey))
      }
      
      winz = session$categ[[category]][c("start","end")][session$categ[[category]][column]==a,]
      if(nrow(winz)>0){
        winz$end[winz$start == winz$end] = winz$end[winz$start == winz$end] +1
        winz = winz[winz$end<end(stream)[1]+1,]
        # if(nrow(winz)==0 || max(winz$end)>end(stream)[1]) stop("C'è un errore grave nei minutaggi, le categorie vanno oltre la durata massima dei segnali")
        
        # winList = apply(winz,1,function(i){as.vector(window(stream,start=i["start"],end=i["end"]))})
        winList = lapply(1:nrow(winz),function(w){as.vector(window(stream,start=winz[w,"start"],end=winz[w,"end"]))})
        merged = na.omit(do.call("c",winList))
        sessions[i,] =  c(i,i, #session #dyadId
                          mean(merged),                                #mean of category
                          sd(merged),                                  #sd of category
                          median(merged),
                          NA,                                          #"aver. marci is meaningless",               #mean(sapply(winList,marciIndex), na.rm=T),
                          mean(winz$end - winz$start, na.rm=T),        #mean duration
                          min(winz$end - winz$start, na.rm=T),         #min duration
                          max(winz$end - winz$start, na.rm=T),         #max duration
                          length(winList),                             #n occurrences
                          quantile(merged, probs=0.25),
                          quantile(merged, probs=0.75),
                          min(merged, na.rm=T),
                          max(merged, na.rm=T)
        )
      } else sessions[i,] = c(i,i,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA)
      sessions[i,1] =paste(attr(session,"sessionId"),attr(session,"dyadId"),sep="_")
      sessions[i,2] =attr(session,"dyadId")
    })
    #print(attr(stream,"name"))
    #print(sessions[i,])
    #}
    if(return.type == "summary"){
      res[[a]]=sessions
    }else{
      res[[a]] = merged
      }
    
  }
  #print(attr(stream,"name"))
  return(res)
}

catRandom <- function(x, signal="RRmean", streamKey="bestCCF", category="PACS", column="CATEGORIA",nRand = 1000) {
  if(!is.DyadExperiment(x) & !is.DyadSession(x))stop("only experiments or sessions get butter")  
  UseMethod("catRandom",x)
}

catRandom.DyadExperiment =function(full, signal="RRmean", streamKey="bestLag", category="PACS", column="CATEGORIA",nRand = 1000){
  
  #se non sai pensare da intelligente, pensa da stupido.
  # column = "CATEGORIA"
  # category="PACS"
  # categKey = "involving"
  # streamKey = "bestLag"
  # signal = "RRmean"
  #  nRand = 1000 
  k=1
  setWin = list()
  for(session in full){
    setWin = c(setWin,catRandom(session,signal,streamKey,category,column,nRand))
  }
  return(setWin)
}



catRandom.DyadSession =function(session, signal="RRmean", streamKey="bestCCF", category="PACS", column="CATEGORIA",nRand = 1000){
  
  #se non sai pensare da intelligente, pensa da stupido.
  # column = "CATEGORIA"
  # category="PACS"
  # categKey = "involving"
  # streamKey = "bestLag"
  # signal = "RRmean"
  #  nRand = 1000 
  #k=1
  
  cat("\r\nAnalyzing session:",paste(attr(session,"sessionId"),attr(session,"dyadId"),sep="_"))
  setWin = list()
  #for(session in full){
  #sess = "03_FR"
  myCat = session$categ$PACS[,c("start", "end", "delta", column)]
  #qui catSummary invece dei valori fissi
  catMin = max(1,min(myCat$delta)) #non voglio avere durate=0
  catMax = max(myCat$delta)
  nMin = max(3,min(sapply(split(myCat,myCat[[column]]),nrow))) #n_cat minime per set (almeno 3)
  nMax = max(sapply(split(myCat,myCat[[column]]),nrow)) #n_cat massime per set 
  
  #toglie i dati reali
  cleanStream = noCatStream(getCCF(session$signals[[signal]], streamKey), session$categ[[category]])
  
  #ora crea i minutaggi random
  myCat$set = NA
  ranCat = myCat[NULL,]
  timeMax = end(cleanStream)[1] - catMax -1
  sym=c("\\","|","/","-"); symvec = rep(1:4,(nRand*nMax/4)+1);  symi =1; cat("\r\n")
  for (i in 1:nRand){
    for(j in 1:round(runif(1,nMin,nMax))){
      rStart = round(runif(1,1,timeMax))
      rDelta = round(runif(1,catMin,catMax))
      ranCat=rbind(ranCat,data.frame(rStart,rStart + rDelta,rDelta,"random",set=i))
      #if(is.na(sym[symvec[symi]]))(print(symi))
      cat("\r",sym[symvec[symi]],"\r ");symi = symi+1
    }
  }
  colnames(ranCat) = c("start", "end", "delta", column,"set")
  
  winz = ranCat[,c("start","end","set")]
  
  ##per test sui dati reali esegui qua:
  #winz = myCat[,c("start","end")]
  #winz$end[winz$end -winz$start ==0] =winz$end[winz$end -winz$start ==0] +1
  ##
  
  
  allWin = apply(winz,1,function(i){as.vector(window(cleanStream, start=i["start"],end=i["end"]))})
  allWinNames = Map(function(mystream,set){lead0(rep(set,length(mystream)))},allWin,winz$set)
  my=data.frame(stream=do.call("c",allWin),set=do.call("c",allWinNames)) #my per ciascun set, corrisponde ad un "merged" di catSummary, ed è ok per calcolare 
  setWin[paste(attr(session,"sessionId"),attr(session,"dyadId"),sep="_")] = list(split(my$stream,my$set)) #setWin è una lista di set, ciascun set corrisponde a 'merged' di una session reale
  #k=k+1
  #}
  return(setWin)
}

#this function sets to NA the values of a stream corresponding to the windows of a dyadCategory object.
#it is used to compare the stream times where there are no categories
noCatStream = function(fullstream,dyadCategory){
  rmCat = dyadCategory[,c("start", "end")] #estrae inizio e fine di ogni finestra categoria
  streamFreq = frequency(fullstream) #salva la frequenza di fullstream
  cleanStream = as.vector(fullstream) #toglie tutti i metadata a fullstream
  #per ciascuna finestra...
  for(i in 1:nrow(rmCat)){
    #volendomi fidare di pastJo, qui imposto a NA tutti i valori compresi in una finestra.
    cleanStream[(rmCat[i,"start"] *streamFreq):(rmCat[i,"end"] *streamFreq)-streamFreq+1] =NA #roba brutta uguale a window
  }
  #recreates the ts
  resTS = ts(na.omit(cleanStream),frequency = streamFreq,start=start(fullstream))
  attr(resTS, "dropped") =   sum(is.na(cleanStream))/length(cleanStream) #how much of the original data was deleted?
  return(resTS)
}


#originally exclude was set to default c("mancante","support","secure base","reflection")

plotCatLag = function(dyadData, signal, category, column, nRand = 0, min_cat=0,
                      exclude=c(), include=c(),sync=c("best","pp")
                      )
{
  
  if(!is.DyadExperiment(dyadData) & !is.DyadSession(dyadData))stop("only experiments or sessions get butter")  
  UseMethod("plotCatLag",dyadData)
}
  


plotCatLag.DyadExperiment = function(EXP, signal, category, column, nRand, min_cat, exclude,sync=c("best","pp"))
{  
    dyads = unique(sapply(EXP,attr,"dyadId"))
    ndyads = length(dyads)
    cat("\r\n-----",ndyads,"dyads found rerouting to ")
    if(ndyads ==1) {
      cat("SEGMENT #1")
      #################################
      ## This chunk manages what to plot
      ## if there is only one dyad
      ## most probably this is a
      ## longitudinal study
      
      
      ###--> the following code aggregates the whole experiment and plots the aggregates by category
      all_sess_lag = lapply(EXP, catSummary, signal, stream="ppLag", category, column, return.type="full")
      all_sess_ccf = lapply(EXP, catSummary, signal, stream="ppCor", category, column, return.type="full")
      all_cat = unlist(lapply(all_sess_lag,names))
      toPlot = names(table(all_cat)[table(all_cat)>min_cat])
      toPlot = toPlot[!toPlot%in%exclude]
      #if(lenght )
      count_vec = as.numeric(table(unlist(sapply(EXP, function(x){x$categ$attref$TIPO}))))
      
      res_lag = res_ccf = vector("list",length(toPlot))
      names(res_lag) = names(res_ccf) = toPlot
      
      #get uncategorized data
      streamClean_lag = do.call("c",lapply(EXP, function(session) {noCatStream(getCCF(session$signals[[signal]], "ppLag"), session$categ[[category]])}))
      streamClean_ccf = do.call("c",lapply(EXP, function(session) {noCatStream(getCCF(session$signals[[signal]], "ppCor"), session$categ[[category]])}))
      random_quantiles_lag =  quantile(streamClean_lag,probs=c(0.25,0.5,0.75))
      random_quantiles_ccf =  quantile(streamClean_ccf,probs=c(0.25,0.5,0.75))
      
      for(categ in toPlot){
        res_lag[[categ]] = c(unlist(lapply(all_sess_lag,function(x){x[[categ]]})))
        res_ccf[[categ]] = c(unlist(lapply(all_sess_ccf,function(x){x[[categ]]})))
        
      }
      ccf_vec = sapply(res_ccf,median)
      
      plotCatLagEngine(res_lag, random_quantiles_ccf, random_quantiles_lag, toPlot,count_vec,ccf_vec,
                       main=paste(attr(EXP[[1]],"dyadId"),":",signal,"lag in",category,"categories"))
      
      ###--> the following code produces a longit plot for each category level
      
      for(categ in toPlot){
        res_lag =  lapply(all_sess_lag,function(x){x[[categ]]})
        count_vec = lapply()
        plotCatLagEngine(res_lag, random_quantiles_ccf, random_quantiles_lag, categ,count_vec,ccf_vec,
                         main=paste(attr(EXP[[1]],"dyadId"),":",signal,"lag in",category,"categories"))
        
      }
      

      ##             end             ##
      #################################
    } else {
      #################################
      ## This chunk manages what to plot
      ## if there are multiple dyads:
      ## the focus can be comparing
      ## different subjects, or just get
      ## the effect of categories
      cat("SEGMENT #2")
      
      ###--> the following code aggregates the whole experiment and plots the aggregates by category
      all_sess_lag = lapply(EXP, catSummary, signal, stream="ppLag", category, column, return.type="full")
      all_sess_ccf = lapply(EXP, catSummary, signal, stream="ppCor", category, column, return.type="full")
      all_cat = unlist(lapply(all_sess_lag,names))
      toPlot = names(table(all_cat)[table(all_cat)>min_cat])
      toPlot = toPlot[!toPlot%in%exclude]
      count_vec = as.numeric(table(unlist(sapply(EXP, function(x){x$categ$attref$TIPO}))))
      
      res_lag = res_ccf = vector("list",length(toPlot))
      names(res_lag) = names(res_ccf) = toPlot
      
      #get uncategorized data
      streamClean_lag = do.call("c",lapply(EXP, function(session) {noCatStream(getCCF(session$signals[[signal]], "ppLag"), session$categ[[category]])}))
      streamClean_ccf = do.call("c",lapply(EXP, function(session) {noCatStream(getCCF(session$signals[[signal]], "ppCor"), session$categ[[category]])}))
      random_quantiles_lag =  quantile(streamClean_lag,probs=c(0.25,0.5,0.75))
      random_quantiles_ccf =  quantile(streamClean_ccf,probs=c(0.25,0.5,0.75))
      
      for(categ in toPlot){
        res_lag[[categ]] = c(unlist(lapply(all_sess_lag,function(x){x[[categ]]})))
        res_ccf[[categ]] = c(unlist(lapply(all_sess_ccf,function(x){x[[categ]]})))
        
      }
      ccf_vec = sapply(res_ccf,median)
      
      plotCatLagEngine(res_lag, random_quantiles_ccf, random_quantiles_lag, toPlot,count_vec,ccf_vec,
                       main=paste(attr(EXP[[1]],"dyadId"),":",signal,"lag in",category,"categories"))
      
      ###--> the following code produces a longit plot for each category level
      
      for(categ in toPlot){
        res_lag =  lapply(all_sess_lag,function(x){x[[categ]]})
        count_vec = lapply()
        plotCatLagEngine(res_lag, random_quantiles_ccf, random_quantiles_lag, categ,count_vec,ccf_vec,
                         main=paste(attr(EXP[[1]],"dyadId"),":",signal,"lag in",category,"categories"))
        
      }
      
      
      
      
      ##             end             ##
      #################################
    }
    #################################
    ## This chunk manages what to plot
    ## in both cases (1 or more dyads)
    cat("\r\n Finally: SEGMENT #3")
    
    
    
    
    
    ##             end             ##
    #################################
    
}




plotCatLag.DyadSession = function(session, signal, category, column, nRand, min_cat, exclude=c(), include=c(),sync=c("best","pp")
                      ) {
  ###debug
  # session = full[[1]]
  # signal = "SC"
  # category= "PACS"
  # column="NUOVACATEGORIA"
  # nRand= 0
  # min_cat=2
  
  sync = match.arg(sync)
  if(sync=="best"){
    syncStream = "bestCCF"
    lagStream = "bestLag"
  } else {
    syncStream = "ppSync"
    lagStream = "ppLag"
  }
  
  if(!is.DyadSession(session))stop("only sessions get butter")  
  capture.output({
    allCat = session$categ[[category]][[column]]
    res = catSummary(session, signal, streamKey =lagStream, category, column, return.type="full")
    res_smr = catSummary(session, signal, streamKey =lagStream, category, column, return.type="summary")
    res_ccf = catSummary(session, signal, streamKey =syncStream, category, column)
  })
    ####dati categorie
    toPlot = names(table(allCat)[table(allCat)>min_cat]) #tiene solo categorie con più di min_cat occorrenze
    if(length(include)>0) toPlot = toPlot[toPlot%in%include] # tiene solo le categorie indicate in include
    if(length(exclude)>0) toPlot = toPlot[!toPlot%in%exclude] # toglie le categorie indicate in exclude
    count_vec = as.vector(sapply(res_smr[toPlot], function(x){x$count})) #crea un vettore con il numero di occorrenze
    ccf_vec = as.vector(unlist(sapply(res_ccf[toPlot], function(x){x$median}))) #estrae le mediane delle ccf

    
  # median (or random) lag
    capture.output({
    if(nRand==0){
      #try to compare the categories with the uncategorized stream, only if there is at least 25% uncategorized stream
      #otherwise compare to the whole stream
      streamClean = noCatStream(getCCF(session$signals[[signal]], lagStream), session$categ[[category]])
      if(attr(streamClean,"dropped")>0.75){
        streamClean = na.omit(getCCF(session$signals[[signal]], lagStream))
      }
      random_quantiles_lag = quantile(streamClean,probs=c(0.25,0.5,0.75))
    } else{
      ranSess = catRandom(session,signal,lagStream,category,column,nRand)
      random_quantiles_lag = quantile(do.call("c",ranSess[[1]]), probs=c(0.25,0.5,0.75))
    }})
    

    #median ccf
    capture.output({
      #try to compare the categories with the uncategorized stream, only if there is at least 25% uncategorized stream
      #otherwise compare to the whole stream
    streamClean = noCatStream(getCCF(session$signals[[signal]], syncStream), session$categ[[category]])
    if(attr(streamClean,"dropped")>0.75){
      streamClean = na.omit(getCCF(session$signals[[signal]], syncStream))
    }
    })
    random_quantiles_ccf = quantile(streamClean, probs=c(0.25,0.5,0.75))
    
    #plot
    plotCatLagEngine(res, random_quantiles_ccf, random_quantiles_lag, toPlot,count_vec,ccf_vec,
                     main=paste(attr(session,"dyadId"),attr(session,"sessionId"),":",signal,"lag in",category,"categories"))
    
}  


plotCatLagEngine = function(res, random_quantiles_ccf, random_quantiles_lag, toPlot,count_vec,ccf_vec,main
) {
  #debug
  cat("\r\nplotCatLagEngine started")
  
  x_at=seq(1,length(toPlot)*2,by=2)
  x_at_lab=seq(1,max(x_at),length.out = length(toPlot)) #general
  rLines = random_quantiles_lag
  medCCF = random_quantiles_ccf[2]
  
  all_col= dual.colors(201,center=medCCF)
  ccf_col = all_col[ as.vector(round(ccf_vec,2)*100 + 101)]
  
  lim=max(sapply(res,function(x) max(abs(x),na.rm=T)))
  
  boxplot(res[toPlot], at=x_at,ylim=c(-lim-1,+lim+1.5),xaxt="n",xlab = "",
          xlim=c(0,max(x_at)+4),ylab="lag (s)",col=ccf_col,
          main=main)
  
  abline(h=0,lty=3)
  abline(h=rLines,col="red",lty=c(2,3,2))
  
  text(c(x_at,max(x_at)+2.5),lim+1.5,labels=c(round(ccf_vec,2),round(medCCF,3) ))
  
  
  rect(max(x_at)+2,   seq(-lim/2,lim/2,by=lim/201),  max(x_at)+3,  seq(-lim/2,lim/2,by=lim/201)+lim/201, col=all_col,border=NA)
  text(max(x_at)+3.3,y=c(-lim/2,0,lim/2,lim/2+2), labels=c("-1","0","+1","SYNC"))
  text(max(x_at)+3.3,y=1.5, labels=paste0(round(medCCF,2)," ± ", round(sd(random_quantiles_ccf),2) ), cex=0.8 )
  segments(x_at_lab+1,-lim-4,x_at_lab+1,lim+4,lwd=0.5)
  axis(side = 1,at = x_at_lab-0.3,tick = F,labels = gsub(" ","\n",toPlot),las=2,cex.axis=0.9)
  axis(side = 1,at = x_at_lab-0.3,tick = F,labels = paste0("\n\n\nn=",count_vec),las=2,cex.axis=1)
  
  # text(x_at+0.4,y=quart3_vec+3,labels = count_vec, cex=0.8)
  # text(x_at+0.4,y=quart3_vec+4, labels =dyadId,cex=0.6)
}  



plotCatLagCompare = function(sessionList, signal="RRmean", category="PACS", column="CATEGORIA", nRand = 0,
                      exclude=c("mancante","support","secure base","reflection")
) {
  if(!is.DyadSession(sessionList))stop("only sessions get butter")  
  
  ##step 1:individuare TUTTI i livelli possibili in tutto l'esperimento
  allCat = factor(as.character(unlist(sapply(sessionList, function(session){as.character(session$categ[[category]][,column])}))))
  allLevels = levels(allCat)
  full2= lapply(full,function(session){
    session$categ[[category]][,column] = factor(as.character(session$categ[[category]][,column]),levels=allLevels)
    session
  })
  #attributes(full2) =attributes(full)
  
  
  streamKey="bestLag"
  
  capture.output({
    allCat = session$categ[[category]][[column]]
    res = catSummary(session, signal, stream="bestLag", category, column, return.type="full")
    res_smr = catSummary(session, signal, stream="bestLag", category, column, return.type="summary")
    res_ccf = catSummary(session, signal, stream="bestCCF", category, column)
  })
  ####dati categorie
  toPlot = names(table(allCat)[table(allCat)>2])
  toPlot = toPlot[!toPlot%in%exclude]
  #print(toPlot)
  #par(oma=c(2,0,0,0))
  
  #print(length(toPlot)*2)
  x_at=seq(1,length(toPlot)*2,by=2)
  x_at_lab=seq(1,max(x_at),length.out = length(toPlot)) #general
  
  count_vec = as.vector(sapply(res_smr[toPlot], function(x){x$count}))
  dyadId = as.vector(unlist(sapply(res_smr[toPlot], function(x){x$dyadId})))
  
  
  ccf_vec = as.vector(unlist(sapply(res_ccf[toPlot], function(x){x$median})))
  
  ##dati random o tutti i dati non categorizzati?
  capture.output({
    if(nRand==0){
      streamClean = noCatStream(getCCF(session$signals[[signal]], "bestLag"), session$categ[[category]])
      rLines = quantile(streamClean,probs=c(0.25,0.5,0.75))
    } else{
      ranSess = catRandom(session,signal,"bestLag",category,column,nRand)
      rLines = quantile(do.call("c",ranSess[[1]]), probs=c(0.25,0.5,0.75))
    }})
  
  #median ccf
  capture.output({
    streamClean = noCatStream(getCCF(session$signals[[signal]], "bestCCF"), session$categ[[category]])
  })
  medCCF = median(streamClean)
  all_col= dual.colors(201,center=medCCF)
  ccf_col = all_col[ as.vector(round(ccf_vec,2)*100 + 101)]
  
  lim=max(sapply(res,function(x) max(abs(x),na.rm=T)))
  
  boxplot(res[toPlot], at=x_at,ylim=c(-lim-1,+lim+1),xaxt="n",xlab = "",
          xlim=c(0,max(x_at)+4),ylab="lag (s)",col=ccf_col,
          main=paste(attr(session,"dyadId"),attr(session,"sessionId"),":",signal,"lag in",category,"categories"))
  
  abline(h=0,lty=3)
  abline(h=rLines,col="red",lty=c(2,3,2))
  
  
  rect(max(x_at)+2,   seq(-lim/2,lim/2,by=lim/201),  max(x_at)+3,  seq(-lim/2,lim/2,by=lim/201)+lim/201, col=all_col,border=NA)
  text(max(x_at)+3.3,y=c(-lim/2,0,lim/2,lim/2+2), labels=c("-1","0","+1","SYNC02"))
  segments(x_at_lab+1,-lim-4,x_at_lab+1,lim+4,lwd=0.5)
  axis(side = 1,at = x_at_lab-0.3,tick = F,labels = gsub(" ","\n",toPlot),las=2,cex.axis=0.9)
  axis(side = 1,at = x_at_lab-0.3,tick = F,labels = paste0("\n\n\nn=",count_vec),las=2,cex.axis=1)
  
  # text(x_at+0.4,y=quart3_vec+3,labels = count_vec, cex=0.8)
  # text(x_at+0.4,y=quart3_vec+4, labels =dyadId,cex=0.6)
}  



#This function splits a stream into category windows and creates boxplots
plotCatStream = function(dyadData,      #a DyadExperiment or DyadSession object
                         signal,        #string indicating the name of the signal
                         stream,        #string indicating the name of the stream
                         category,      #string of the name of the DyadCategory object
                         column,        #string of the name of a factor column in the DyadCategory
                         nRand = 0,     #if to compare with random iteration, how many? If 0, either the uncategorized stream or the whole stream will be used for intervals
                         least.cat.n=0, #minimum number of instances for a factor level to be included
                         exclude=c()    #explicitly exclude factor levels by name
)
{
  
  if(!is.DyadExperiment(dyadData) & !is.DyadSession(dyadData))stop("only experiments or sessions get butter")  
  UseMethod("plotCatStream",dyadData)
}

plotCatStream.DyadExperiment = function(EXP, signal, stream, category, column, nRand, least.cat.n, exclude)
{  
  dyads = unique(sapply(EXP,attr,"dyadId"))
  ndyads = length(dyads)
  cat("\r\n-----",ndyads,"dyads found rerouting to ")
  if(ndyads ==1) {
    cat("Longitudinal design")
    #################################
    ## This chunk manages what to plot
    ## if there is only one dyad
    ## most probably this is a
    ## longitudinal study
    
    
    ###--> the following code aggregates the whole experiment and plots the aggregates by category
    ###
    
    #split each session by category with catSummary
    all_sess_stream = lapply(EXP, catSummary, signal, stream=stream, category, column, return.type="full")
    #get real category names
    all_cat = unlist(lapply(all_sess_stream,names))
    #apply least.cat.n and exclude criteria
    toPlot = names(table(all_cat)[table(all_cat)>least.cat.n])
    toPlot = toPlot[!toPlot%in%exclude]
    ncat = length(toPlot)
    #if(lenght )
    count_vec = table(unlist(sapply(EXP, function(x){x$categ[[category]][[column]]})))
    
    #initialize empty result vector
    res  = vector("list",length(toPlot))
    names(res) = toPlot
    
    #clean: get data not corresponding to any category
    streamCleanList =lapply(EXP, function(session) {noCatStream(getCCF(session$signals[[signal]], stream), session$categ[[category]])})
    dropped = sapply(streamCleanList, attr, "dropped") #check how much signal is lost through cleaning
    if (mean(dropped, na.rm=T)>0.80 | max(dropped)>0.80){
      warning("The categories extended over 80% of signal. The whole experiment will be used for comparison")
      streamCleanList = do.call("c",all_sess_stream)
    }
    streamClean = do.call("c",streamCleanList)
    random_quantiles =  quantile(streamClean,probs=c(0.25,0.5,0.75))
    
    for(categ in toPlot){
      res[[categ]] = c(unlist(lapply(all_sess_stream,function(x){x[[categ]]})))
    }
    
    # ccf_vec = sapply(res_ccf,median)
    # par(mfrow=c(1,1))
    plotCatStreamEngine(res, random_quantiles, toPlot,count_vec, ylab = stream,
                        main=paste(attr(EXP[[1]],"dyadId"),":",signal,stream,"in",category,column,"categories"))
    
    ###--> the following code produces a longit plot for each category level
    
    # if(ncat>1){
    #   if(ncat==2) par(mfrow=c(2,1))
    #   else par(mfrow=c(2,2))
    # }
    for(categ in toPlot){
      resx =  lapply(all_sess_stream,function(x){x[[categ]]})
      boxplot(resx, main=categ,xaxt="n",ylab=stream)
      axis(side = 1,at = seq_along(all_sess_stream),labels = sapply(EXP, attr, "sessionId"),las=2,cex.axis=0.9)
      abline(h=0,lty=3)
      abline(h=random_quantiles,col="red",lty=c(2,3,2))
      
      # plotCatLagEngine(res_lag, random_quantiles_ccf, random_quantiles_lag, categ,count_vec,ccf_vec,
      #                  main=paste(attr(EXP[[1]],"dyadId"),":",signal,"lag in",category,"categories"))
      
    }
    # par(mfrow=c(1,1))
    
    
    ##             end             ##
    #################################
  } else {
    #################################
    ## This chunk manages what to plot
    ## if there are multiple dyads:
    ## the focus can be comparing
    ## different subjects, or just get
    ## the effect of categories
    cat("Comparison of cases - not written yet :-(")
    
    # ###--> the following code aggregates the whole experiment and plots the aggregates by category
    # all_sess_lag = lapply(EXP, catSummary, signal, stream="bestLag", category, column, return.type="full")
    # all_sess_ccf = lapply(EXP, catSummary, signal, stream="bestCCF", category, column, return.type="full")
    # all_cat = unlist(lapply(all_sess_lag,names))
    # toPlot = names(table(all_cat)[table(all_cat)>least.cat.n])
    # toPlot = toPlot[!toPlot%in%exclude]
    # count_vec = as.numeric(table(unlist(sapply(EXP, function(x){x$categ$attref$TIPO}))))
    # 
    # res_lag = res_ccf = vector("list",length(toPlot))
    # names(res_lag) = names(res_ccf) = toPlot
    # 
    # #get uncategorized data
    # streamClean_lag = do.call("c",lapply(EXP, function(session) {noCatStream(getCCF(session$signals[[signal]], "bestLag"), session$categ[[category]])}))
    # streamClean_ccf = do.call("c",lapply(EXP, function(session) {noCatStream(getCCF(session$signals[[signal]], "bestCCF"), session$categ[[category]])}))
    # random_quantiles_lag =  quantile(streamClean_lag,probs=c(0.25,0.5,0.75))
    # random_quantiles_ccf =  quantile(streamClean_ccf,probs=c(0.25,0.5,0.75))
    # 
    # for(categ in toPlot){
    #   res_lag[[categ]] = c(unlist(lapply(all_sess_lag,function(x){x[[categ]]})))
    #   res_ccf[[categ]] = c(unlist(lapply(all_sess_ccf,function(x){x[[categ]]})))
    #   
    # }
    # ccf_vec = sapply(res_ccf,median)
    # 
    # plotCatLagEngine(res_lag, random_quantiles_ccf, random_quantiles_lag, toPlot,count_vec,ccf_vec,
    #                  main=paste(attr(EXP[[1]],"dyadId"),":",signal,"lag in",category,"categories"))
    # 
    # ###--> the following code produces a longit plot for each category level
    # 
    # for(categ in toPlot){
    #   res_lag =  lapply(all_sess_lag,function(x){x[[categ]]})
    #   count_vec = lapply()
    #   plotCatLagEngine(res_lag, random_quantiles_ccf, random_quantiles_lag, categ,count_vec,ccf_vec,
    #                    main=paste(attr(EXP[[1]],"dyadId"),":",signal,"lag in",category,"categories"))
    #   
    # }
    # 
    # 
    # 
    
    ##             end             ##
    #################################
  }
  #################################
  ## This chunk manages what to plot
  ## in both cases (1 or more dyads)
  cat("\r\n Finally: SEGMENT #3")
  
  
  
  
  
  ##             end             ##
  #################################
  
}


# 
# 
# plotCatStream.DyadSession = function(session, signal, category, column, nRand, least.cat.n, exclude
# ) {
#   if(!is.DyadSession(session))stop("only sessions get butter")  
#   capture.output({
#     allCat = session$categ[[category]][[column]]
#     res = catSummary(session, signal, stream="bestLag", category, column, return.type="full")
#     res_smr = catSummary(session, signal, stream="bestLag", category, column, return.type="summary")
#     res_ccf = catSummary(session, signal, stream="bestCCF", category, column)
#   })
#   ####dati categorie
#   toPlot = names(table(allCat)[table(allCat)>least.cat.n]) #tiene solo categorie con più di least.cat.n occorrenze
#   toPlot = toPlot[!toPlot%in%exclude] # toglie le categorie indicate in exclude
#   count_vec = as.vector(sapply(res_smr[toPlot], function(x){x$count})) #crea un vettore con il numero di occorrenze
#   ccf_vec = as.vector(unlist(sapply(res_ccf[toPlot], function(x){x$median}))) #estrae le mediane delle ccf
#   
#   
#   # median (or random) lag
#   capture.output({
#     if(nRand==0){
#       #try to compare the categories with the uncategorized stream, only if there is at least 25% uncategorized stream
#       #otherwise compare to the whole stream
#       streamClean = noCatStream(getCCF(session$signals[[signal]], "bestLag"), session$categ[[category]])
#       if(attr(streamClean,"dropped")>0.75){
#         streamClean = na.omit(getCCF(session$signals[[signal]], "bestLag"))
#       }
#       random_quantiles_lag = quantile(streamClean,probs=c(0.25,0.5,0.75))
#     } else{
#       ranSess = catRandom(session,signal,"bestLag",category,column,nRand)
#       random_quantiles_lag = quantile(do.call("c",ranSess[[1]]), probs=c(0.25,0.5,0.75))
#     }})
#   
#   
#   #median ccf
#   capture.output({
#     #try to compare the categories with the uncategorized stream, only if there is at least 25% uncategorized stream
#     #otherwise compare to the whole stream
#     streamClean = noCatStream(getCCF(session$signals[[signal]], "bestCCF"), session$categ[[category]])
#     if(attr(streamClean,"dropped")>0.75){
#       streamClean = na.omit(getCCF(session$signals[[signal]], "bestCCF"))
#     }
#   })
#   random_quantiles_ccf = quantile(streamClean, probs=c(0.25,0.5,0.75))
#   
#   #plot
#   plotCatLagEngine(res, random_quantiles_ccf, random_quantiles_lag, toPlot,count_vec,ccf_vec,
#                    main=paste(attr(session,"dyadId"),attr(session,"sessionId"),":",signal,"lag in",category,"categories"))
#   
# }  


plotCatStreamEngine = function(res, random_quantiles, toPlot, count_vec, main, ylab
) {
  #debug
  cat("\r\nplotCatLagEngine started")
  print(toPlot)
  x_at=seq(1,length(toPlot)*2,by=2)-0.5
  x_at_lab=seq(1,max(x_at)+0.5,length.out = length(toPlot))-0.5 #general
  rLines = random_quantiles
  
  # all_col= dual.colors(201,center=medCCF)
  # ccf_col = all_col[ as.vector(round(ccf_vec,2)*100 + 101)]
  
  lim1=max(unlist(res,recursive = T))
  lim2=min(unlist(res,recursive = T))
  
  boxplot(res[toPlot], at=x_at,ylim=c(lim2,lim1),xaxt="n",xlab = "",
          xlim=c(0,max(x_at)+0.5),ylab=ylab,
          main=main)
  
  abline(h=0,lty=3)
  abline(h=rLines,col="red",lty=c(2,3,2))
  
  segments(x_at_lab+1,-lim1-4,x_at_lab+1,lim2+4,lwd=0.5)
  axis(side = 1,at = x_at_lab-0.3,tick = F,labels = gsub(" ","\n",toPlot),las=2,cex.axis=0.9)
  axis(side = 1,at = x_at_lab-0.3,tick = F,labels = paste0("\n\n\nn=",count_vec),las=2,cex.axis=1)
  
}  





### da qua sotto roba vecchia  
#   
#   i=0
#   for(a in res[toPlot]){
#     i=i+1
#     
#     #png(paste0(plotPath,"\\",paste(signal,streamKey,category,toPlot[i],collapse="_",sep="_"),".png"), width = 800, height = 800, units = "px" )
#     par(mfrow=c(2,2),oma=c(0,0,2,0))
#     
#     myYlim = c(-1,1)
#     cols = c("dodgerblue4", "deeppink3")
#     if("ylim"%in%names(dots)) {
#       print("TRUE")
#       meanylim = dots$ylim
#       #dots$ylim = NULL
#     } else  meanylim = myYlim
#     print(meanylim)
#     barplot(a[,2], main="mean sync", ylim = meanylim,names.arg = a[,1],col = cols)
#     barplot(a[,3],main="sd sync", ylim = myYlim,names.arg = a[,1],col = cols)
#     barplot(a[,4],main="marci index", ylim = myYlim*3,names.arg = a[,1],col = cols)
#     barplot(a[,5],main="duration (s)",names.arg = a[,1],col = cols,ylim =c(0,60))
#     title(main=paste0(toPlot[i]," (",paste(paste(a[,1],a[,6],sep=": "),collapse="; "),")"),outer=T)
#     title(main=paste0("signal: ", signal,"; stream: ",streamKey),outer=T,cex.main=0.8, line=-1)
#     
#     #dev.off()
#   }
#   return(res)
# }
# 
# 
# #catSummary(full, signal="SC", stream="bestCCF", category="PACS", column="CATEGORIA")
# 
# 
# 
# 
# 
# 
# 
# 
# #this function extracts a column from ccfmat of a given signal
# #the data is extracted as a DyadStream object with appropriate frequency
# getCCF = function(signal, lag, col="red", lty=2,lwd=2){
#   if(!is.DyadSignal(signal)) stop("Only objects of class DyadSignal can be processed by this function")
#   if(!grepl("lag",lag, ignore.case = T) & grepl("best",lag, ignore.case = T)) {
#     lag = "bestCCF"
#   }  else if (grepl("lag",lag, ignore.case = T) & grepl("best",lag, ignore.case = T)) {
#     lag = "bestLag"
#   } else if(suppressWarnings(!is.na(as.numeric(lag)))){
#     lag = paste0("lag",as.character(lag))
#   }
#   if (! lag %in% colnames(signal$ccf$ccfmat)) stop("Specified lag value '",lag,"' was not found. Available lags: ", paste(colnames(signal$ccf$ccfmat), collapse=" | "))
#   cat0("\r\n CCF column '",lag,"' was selected")
#   if(!signal$ccf$settings$interpolated) warning("CCF is not interpolated. Sample rate setting might be wrong, or manifest other unexpected result",call.=F)
#   
#   stream = DyadStream(ts(signal$ccf$ccfmat[as.character(lag)], frequency = signal$ccf$sampRate),name = paste0("CCF - ",lag), col = col, lty = lty, lwd = lwd)
#   CCFStream(stream, 
#             lagSec       = signal$ccf$settings$lagSec,
#             incSec       = signal$ccf$settings$incSec,
#             winSec       = signal$ccf$settings$winSec,
#             accelSec     = signal$ccf$settings$accelSec,
#             weight       = signal$ccf$settings$weight,
#             interpolated = signal$ccf$settings$interpolated)
# }
