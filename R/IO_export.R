
############################################################################################################
############################################################################################################
############################################################################################################
## The final exporter
# Note: use selectSignals() to reduce the number of columns



#' Title
#'
#' @param dirPath 
#' @param experiment 
#' @param signals 
#' @param CCF 
#' @param onlyBest 
#' @param original 
#'
#' @return
#' @export
#'
#' @examples


megaExpExport = function(dirPath, experiment, signals="all", CCF=T, onlyBest=F, original=T){
  warning("This function is VEEEEERY old in development queue, and was not tested since 2016")
  ###DEBUG
  #   dirPath = "__RLab\\demo_data\\expExp"
  #   experiment = qwe2
  #   signals = "all"
  #   CCF =T
  #   onlyBest = F
  #   original = T
  ####
  # rm(signal,session,iSession,original,onlyBest,CCF,signals,experiment,dirPath)
  ####
  
  if(!is(experiment,"DyadExperiment")) stop("Only objects of class DyadExperiment can be processed by this function")
  if(!file.exists(dirPath)){stop("Selected directory does not exists")}
  if(sum(CCF, onlyBest,original) == 0) stop("at least one export parameter must be set TRUE")
  nSessions  = length(experiment)
  Map(function(session,iSession){
    ##debug
    # session = experiment$sessions[[1]]
    # iSession = 1
    ##
    prog(iSession,nSessions)
    if(signals =="all") signals = names(session)
    signalListSR = as.vector(unlist(sapply(session[signals],function(x){c(x$sampRate,x$ccf$sampRate)})))
    #print(str(signalListSR))
    if(!all(signalListSR == signalListSR[1])){
      stop("Not all the output signals have the same sampling rate, please use interpolation!")
    }
    sampRate = signalListSR[1]
    signalList = lapply(session[signals], function(signal){
      ##debug
      # signal = session[[8]]
      ##
      #print(name(signal))
      #rm(signal,my.ccf,my.orig,uberDS)
      if(is(signal,"DyadSignal")){
        if(CCF && !is.null(signal$ccf)){
          if(onlyBest){
            my.ccf = signal$ccf$ccfmat[,c("lag0","bestCCF","bestLag")]
            colnames(my.ccf) = paste(name(signal),c("lag0","bestCCF","bestLag"),sep="_")
          } else {
            my.ccf = signal$ccf$ccfmat
            colnames(my.ccf) = paste(name(signal),colnames(signal$ccf$ccfmat),sep="_")
          }
        } else {
          my.ccf= NULL
        }
        if(original){
          my.orig = data.frame(cbind( signal$s1, signal$s2 ))
          colnames(my.orig) = c(paste(name(signal),attr(signal$s1,"name"),sep="_"), paste(name(signal),attr(signal$s2,"name"),sep="_"))
        }
        unequalCbind(my.orig,my.ccf)
        
      } else if(is(signal,"DyadCategory")){
        my.cat = data.frame(signal$value)
        colnames(my.cat) = name(signal)
        my.cat
      }
    })
    uberDS = do.call("unequalCbind",signalList)
    write.table(uberDS,paste0(dirPath,"\\",session$sessionId,"_",session$dyadId,"_",sampRate,"Hz.csv"),sep=";",fileEncoding = "UTF-8",row.names = F)
  },experiment, seq_along(experiment))
  return(NULL)
}
