#' This library provides tools to overlay signals on video
#'
#' to allow for multiple signals, everything must be a list of length n
#' so if I have e.g. SC, HR and MEA, I must provide 3 dyadsignals, 3 y positions, etc
#' 
#' Directory organization
#' Each of these endeavors require:
#' - json file of the original acquisition
#' - csv version for DyadSync import
#' - video of the session in a suitable format
#' - folders with all the frames
#' - output video file
#' THUS it would make sense to have a folder for each session instead of a folder
#' for each 'type' of data

#' Draws video frames showing dyadic signals
#'
#' @param session 
#' @param signals 
#' @param sync 
#' @param WIN 
#' @param outputDir 
#' @param w 
#' @param h 
#' @param frameRate 
#' @param layout 
#' @param col_sync 
#' @param ffmpeg override the default ffmpeg instruction. Use {IN} to mark the dynamic
#' original video location. Not yet implemented
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
drawFrames = function(session, signals, sync, WIN, inputDir,
                      outputDir, idOrder=c(F,F,T), idSep="_", 
                      w, h, frameRate,
                      layout = c(1), col_sync,
                      ffmpeg.query, ffmpeg.path, multi.core = TRUE, 
                      speed=1, keepAudio = TRUE, ... ){
  ############
  # session = lr[[1]]
  # signals = "SC"
  # sync="amico1"
  # WIN=30
  # inputDir = input.video.path
  # outputDir = file.path("C:","videoConversion","video_sync","vetere")
  # idSep = "_"
  # idOrder = c(F,F,T)
  # w = 640
  # h = 480
  # frameRate = 10
  # ffmpeg.path=file.path("C:","videoConversion","ffmpeg")
  # multi.core = 15
  # speed=1
  # keepAudio = TRUE
  # # # session$HR = session$SC
  # # layout = 1
  # # w = 640#1920#1280#960
  # # h = 480#1080#720#540
  # # WIN = 30 #how much signal is shown at any given point
  # # frameRate = 10
  # # col_sync= c("#f70c3f","#c16107", "#878787","#878787","#878787", "#86a817", "#23c647")
  # # sync = "amico1"
  # # signals = c("SC")#,"HR")
  # # outputDir = file.path("C:","videoConversion","video_sync","vetere")
  # # multi.core= TRUE
  # # ffmpeg.path = file.path("C:","videoConversion","ffmpeg")
  # # speed=1
  # stop("debug")
  #############
  # if(missing(ffmpeg)){
  #   
  # } else {
  #   
  # }
  
  # Setup
  win = round(WIN*frameRate)
  win2 = trunc(win/2)
  iNameF  = UID(session)
  
  #paths THIS PROBABLY BE CHANGED
  sd = file.path(outputDir, iNameF) #session directory
  fd = file.path(sd,"frames")
  #a better logic probably is: IF folder is there, skip frame generation and only do ffmpeg
  
  if(!dir.exists(sd)){
    dir.create(fd, showWarnings = T, recursive=T)
    # stop(paste(sd, "was just created. Please put the video file inside it."))
  } 
  
  if(dir.exists(fd)){
    unlink(fd, recursive = TRUE)
  }
  dir.create(fd, showWarnings = T, recursive=T)
  

  #check video files
  inputVideoName = list.files(inputDir, pattern = ".mp4",full.names = FALSE)
  sessions = sapply(strsplit(inputVideoName,paste0("[\\.,",idSep,"]")),\(x)x[idOrder] ) 
  sessions = sapply(sessions, \(x)x[1])
  sessions = as.numeric(sessions)
  if(any(is.na(sessions))) stop ("sessions could not be extracted from", paste(inputVideoName[which(is.na(sessions))]))
  ss = which(sessions == sessionId(session))
  # if(length(inputVideoName)>1) stop("Only one video file must be present in the folder")
  inputVideo =  list.files(inputDir, pattern = ".mp4",full.names = TRUE)[ss]

  outputVideo = file.path(sd, paste0(iNameF,"_",sync,".mp4"))
  

  # check ffmpeg installation folder
  bin = file.path(ffmpeg.path,"bin")
  if(!file.exists(file.path(bin,"ffmpeg.exe")) ||
     !file.exists(file.path(bin,"ffprobe.exe"))) stop("up to date ffmpeg must be installed in", ffmpeg.path, "with .exe in the bin subdirectory")
  probe = paste0(file.path(bin,"ffprobe.exe -show_streams"),' "',inputVideo,)
  ffmpeg.exe = file.path(bin,"ffmpeg")
  
  ################
  ##Check video properties
  
  #extrapolate average video frameRate
  x = system("cmd.exe",input=probe, intern=TRUE, show=FALSE)
  avg = sapply(x, grepl, pattern="avg_frame_rate" )
  avg = which(avg)
  # Now I am assuming that the first stream is always the video
  avg = x[avg[1]]
  avg = strsplit(avg, "=")[[1]][2]
  avg = strsplit(avg, "/")[[1]]
  avg = as.numeric(avg[1]) / as.numeric(avg[2])
  
  vidur =  sapply(x, grepl, pattern="duration=" , fixed =TRUE)
  vidur = which(vidur)
  vidur = x[vidur[1]]
  vidur = strsplit(vidur, "=")[[1]][2]
  vidur = as.numeric(vidur)
  if(is.na(vidur)) stop("Could not determine video duration :-(")
  

  # Setup color scale for synchrony
  if(missing(col_sync)) col_sync= c("#f70c3f","#c16107", "#878787","#878787","#878787", "#86a817", "#23c647")
  colfunc <- grDevices::colorRampPalette(colors=col_sync, bias=1)
  legendSteps =20
  colz = colfunc(legendSteps)
  bins = seq(1/legendSteps,1, by=1/legendSteps)


  

  # if mp4 video is there, skip entirely
  # if(dir.exists(sd)){
  #   #this directory exists. Do something special
  #   ui = readline("Directory already exists. Do you want to overwrite the data in it? Type yes or no: ")
  #   ui = substr(ui, 1, 1)
  #   while(!ui %in% c("y","n")){
  #     ui = readline("Only y or n are acceptable answers. Retry: ")
  #     ui = substr(ui, 1, 1)
  #   }
  #   if(ui == "n") {
  #     return(paste("Directory", iNameF, "already exists, skipping.")) 
  #     } else {
  #     dir.create(sd, showWarnings = T, recursive=T)
  #     dir.create(fd, showWarnings = T, recursive=T)
  #     
  #   }
  # 
  # } 
  

  

  if(!all(signals %in% names(session))) stop("Specified signals not found in session ", iNameF)
  #restrict session to used signals
  s = session[signals] 
  #resample to video frameRate
  s = signalFilter(s, DyadSync::resample, newSampRate = frameRate)
  
  
  start_values = sapply(s, start)
  vstart = min(start_values)
  end_values = sapply(s, end)
  vend = max(end_values)
  
  
  ## NB: il video inizia sempre a 00:00
  #' quindi il tempo tra 00:00 e start(signal) va riempito di frame bianchi
  #' se start(signal) è == 00:00 bisogna aggiungere WIN/2 sample perché il primo sample
  #' deve comparire al centro della finestra
  
  frame_start = min(0, vstart - WIN/2)
  frame_end   = max(vidur, vend + WIN/2)
  
  # npadframes = round(vstart * frameRate)
  # 
  # 
  # #the signal frames correspond to the resampled data
  # lenz = sapply(s, \(x) length(x$s1))
  # nframes = max(lenz)
  # # TODO
  # #to allow roll-in the first window must have the first sample as the last value
  # #so you need to add win-1 sample to the start (end end) of the signal
  # #empty padding frames (before start) + lenz signal frames + 1 full window to allow roll out
  # 
  # # allframes = npadframes + nframes 
  # allframes = npadframes + nframes + 2*(win)
  
  ## VIDE

  ##Add half window at the start and end of each signal, and uniform the start values
  signal = signals[1]
  sx = s  
  
  nsignals = length(signals)
  signal_height = 0.3
  low_b = cumsum(rep(signal_height,nsignals)) -signal_height + 0.05
  high_b =  cumsum(rep(signal_height,nsignals)) + 0.05
  for(ii in 1:length(sx)){
    log = capture.output({
      sx[[ii]]$s1 = rescaleByWin(sx[[ii]]$s1,WIN=60*5, newMin=low_b[ii], newMax=high_b[ii], cores=1)
      sx[[ii]]$s2 = rescaleByWin(sx[[ii]]$s2,WIN=60*5, newMin=low_b[ii], newMax=high_b[ii], cores=1)
    })

  }
  # window(s$SC$s1,  start = 0, end = 2)
  # window(sx$SC$s1,  start = 0, end = 2)
  # window(s$SC$s1,  start = 4934.1, end = 4964.1)
  # window(sx$SC$s1,  start = 4934.1, end = 4964.1)
  # window(s$SC$s1,  start = 4800.8, end = 4964.1)
  # window(sx$SC$s1,  start = 4800.8, end = 4964.1)
  # 
  
  
  sx[[signal]]$s1           =  window(sx[[signal]]$s1,           start = frame_start, end = frame_end)
  sx[[signal]]$s2           =  window(sx[[signal]]$s2,           start = frame_start, end = frame_end)
  sx[[signal]][[sync]]$sync =  window(sx[[signal]][[sync]]$sync, start = frame_start, end = frame_end)
  sx[[signal]][[sync]]$lag  =  window(sx[[signal]][[sync]]$lag,  start = frame_start, end = frame_end)

  
  ###test
 
  
  #swinz describes the windows to be plotted
  # swinz = nwin(rats(-(win-1):(nframes+win-1),start=vstart), WIN=WIN, INC = 1/frameRate, flex=FALSE, SR = frameRate,return = "all")
  # swinz = nwin(rats(1:(npadframes + nframes+win-1),start=vstart), WIN=WIN, INC = 1/frameRate, flex=FALSE, SR = frameRate,return = "all")
  swinz = nwin(sx[[signal]]$s1, WIN=WIN, INC = 1/frameRate, flex=FALSE, SR = frameRate,return = "all")
  showtimes = timeMaster(trunc(swinz$mid_t), out = "hour")
  swinz$show_t = showtimes
  
  nframes = nrow(swinz)

  # head(swinz,200)
  # tail(swinz,10)

 
  #rescale synchrony data
  # @HACK here you are using just one synchrony, while allowing for multiple signals
  
  #lag is stored in samples in amico1 so:
  lagResolution = frequency(session[[1]]$s1) 
  sx[[signal]][[sync]]$sync = resample(sx[[signal]][[sync]]$sync, frameRate)
  sx[[signal]][[sync]]$lag  = resample(sx[[signal]][[sync]]$lag, frameRate)
  
  sy = rangeRescale(sx[[signal]][[sync]]$sync,0,1,-1,1)
  LAGs = sx[[signal]][[sync]]$lag / lagResolution
  sy[is.na(sy)] = 0
  max_lag = attr(sx[[signal]][[sync]], "lagSec") 
  
  

  
  # plot(rangeRescale(sx$SC$s1[9000:9604], low_b[1],high_b[1]) )
  

  # lines(sx$SC$s1[9000:9604],col=2)
  


  ####### SETUP PARALLELIZATION
  # progresbar
  pb <- progress::progress_bar$new(
    format = "Drawing frames::percent [:bar] :elapsed | ETA: :eta",
    total = nframes,    # number of iterations
    width = 60, 
    show_after=0 #show immediately
  )
  progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  } 
  opts <- list(progress = progress)
  
  #n of cores
  if(multi.core){
    cores=parallel::detectCores()[1]
    if(is.logical(multi.core)){
      #nothing to do
    } else if(is.numeric(multi.core)){
      cores = min(cores,multi.core )
    }
  } else cores = 1
  
  #initialize
  cl <- parallel::makeCluster(cores) #not to overload your computer
  doSNOW::registerDoSNOW(cl)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  pb$tick(0)
  j=1
  #do the job
  resList <- foreach::foreach(
    j=1:nframes, .options.snow = opts, .errorhandling='stop'#.errorhandling='pass'
  )  %dopar% {
  # for(j in 1:nframes){ ##ONLY FOR DEBUG, UNCOMMENT THIS TO RUN AS SIMPLE LOOP
    #draw blank frames
    
    # ## j iterates on all frames (pad + signal)
    # ## k iterates only on signals
    # if(j <= npadframes) {
    # 
    #     png( file.path(fd,paste0("frame_",DyadSync::lead0(j,7),".png")),width=w, height = h, units = "px",type = "cairo")
    #     par(bg=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
    #     plot(-100, xlim=c(0,1),ylim=c(0,1), axes = F, frame.plot = F, xlab = "",ylab="")
    #     points(0.76, 0.953, cex = 3, pch = 16, col="#ff1111")
    #     text(0.95, 0.95, showtimes[j],cex=2,pos = 2)
    #     dev.off()
    # 
    # } else {
    #   # k = npadframes + j
    #   k = j - npadframes
    
    
    ##ERRORS: frame_0048123.png
      
      png( file.path(fd,paste0("frame_",DyadSync::lead0(j,7),".png")),width=w, height = h, units = "px",type = "cairo")
      par(bg=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
      plot(-100, xlim=c(0,1),ylim=c(0,1), axes = F, frame.plot = F, xlab = "",ylab="", xaxs = "i", yaxs = "i")
      points(0.76, 0.953, cex = 3, pch = 16, col="#ff1111")
      text(0.95, 0.95, showtimes[j],cex=2,pos = 2)
      # text(0.95,0.95, pos=2, labels = "")
      
      for(ii  in 1:nsignals){
        si1 = window(sx[[ii]]$s1, start = swinz$start_t[j], end = swinz$end_t[j])
        # test = window(s$SC$s1, start = swinz$start_t[j], end = swinz$end_t[j])
        si2 = window(sx[[ii]]$s2, start = swinz$start_t[j], end = swinz$end_t[j])
        syi = window(sy        , start = swinz$start_t[j], end = swinz$end_t[j])
        lagi= window(LAGs      , start = swinz$start_t[j], end = swinz$end_t[j])
        
        xvals = seq(0.05,0.95,length.out = win)
        lines(xvals, si1$y, col=2,lwd=5)
        lines(xvals, si2$y, col=4,lwd=5)
        
        #il valore di lag in secondi
        wlag = DyadSync::at(LAGs, time = swinz$mid_t[j])$y
        if(!is.na(wlag)){
        #il valore di s2 corrispondente
        lagMatch =  DyadSync::at(si2, time = swinz$mid_t[j] + wlag)
        
        #la posizione sulle x riscalate
        lagMatch_x = xvals[which(si2$x == lagMatch$x)]

          #       # if(swinz$mid_s[j]>0 && swinz$mid_s[j]<=length(si1)){
          # #il valore di lag, in secondi, al tempo mid_s. VERIFIED
          # wlag = LAGs[swinz$mid_s[j]]$y 
          # # lagi[win2+1]
          # 
          # #to verify you can sum the s2 time and the lag sequence
          # #cropping is needed due to rescaleByWin losing the tail
          # # s2laggedtime = sx[[ii]]$s1$x + LAGs$y[1:length(sx[[ii]]$s1$x)]
          # # s2laggedtime[swinz$mid_s[j]]
          # 
          # #s2 al tempo mid_t + wlag
          # lagMatch = window(sx[[ii]]$s2, start = swinz$mid_t[j]+wlag, duration=10)[1]
          # lagMatch_x = xvals[which(lagi$x == lagMatch$x)]
          # if(length(lagMatch_x)>0){
            lagMatch_y = lagMatch$y
            
            if(!is.na(lagMatch_y)){
              mid_s = win2 + 1 #the first sample of the second half-window
              xmid = xvals[mid_s]
              segments(xmid,si1[win2+1],lagMatch_x,lagMatch_y,lwd=8,col="white")
              segments(xmid,si1[win2+1],lagMatch_x,lagMatch_y,lwd=5,col="grey10")
              
              points(xmid,si1[mid_s],pch=18, col="white",cex=4)
              points(xmid,si1[mid_s],pch=18, col="grey10",cex=3)
              
              
              
              #WIN è mappato su 0.90 unitless width units (UWU)
              #quindi per trasformare il lag da secondi a UWU devi
              # UWU_s = 0.9/WIN #1 secondo è uguale a UWU_s
              points(lagMatch_x,lagMatch_y,pch=18, col="white",cex=4)
              points(lagMatch_x,lagMatch_y,pch=18, col="grey10",cex=3)
              
            }
            


            
            
            #sync
            sync_j = sy[swinz$mid_s[j]]
            # sync_j = at(sy, time = swinz$mid_t[j])
            
            polygon(x = c(0.95,1.0,1.0,0.95),y =c(0,0,sync_j+0.01,sync_j+0.01 ), col=colz[sum(sync_j>=bins)+1],border = NA)
            
          # }
          

        }
       
        
        

        text(0,0.95,paste0(names(sx)[1]," - ", s1Name(sx[[1]])),cex=3,col=2)
        text(0,0.90,paste0(names(sx)[1]," - ", s2Name(sx[[1]])),cex=3,col=4)
        
        abline(h = low_b)
        abline(h=high_b)
      }
      dev.off()
    # }
  }
  parallel::stopCluster(cl); graphics.off();beepr::beep()
  



  # setwd("C:/videoConversion")
  
  
  
  complexFilterVideo = paste(
  # "[0:a]asetpts=",1/speed,"*PTS[resa];",
  "[0:v]scale=2*trunc(iw*sar/2):ih,setsar=1[0v_1];",
  # "[0:a]atempo=2.0[resa]",
  "[1:v]format=yuva420p,colorchannelmixer=aa=0.8[ckout];",
  # "[1:v]format=argb,geq=r='r(X,Y)':a='0.5*alpha(X,Y)'[ckout];",
  
  "[0v_1][ckout]overlay[resv];",
  paste0("[resv]setpts=",1/speed,"*PTS[resv];"),
  "[0:a]atempo=",speed,"[resa];",
  
  "")

  
  
  # if(keepAudio){
  #   ffmpeg.audio = "-map 0:a -c:a copy" 
  # } else {
  #   ffmpeg.audio = ""
  # }
  # 
  
    
  
  ffmpeg = paste0(ffmpeg.exe, ' ',
                  '-y ',
                  '-thread_queue_size 1024 ',
                  '-i "',inputVideo,'" -framerate ',frameRate,' ',
                  '-i "',fd,'\\frame_%07d.png" -r ',avg,' ',#-pix_fmt yuva420p ',
                  # '-f lavfi -t 0.1 -i anullsrc',
                  '-filter_complex "',complexFilterVideo,'"  ',
                  # '-filter_complex "',complexFilterAudio,'"  ',
                  '-map "[resv]" -c:v libx264 ',
                  '-map "[resa]" ',
                  '"',outputVideo,'"')
  
  # ffmpeg = paste0('cd C:\\videoConversion && ffmpeg -y -thread_queue_size 1024 -i "video_sync\\vetere\\Vetere Seduta ',iNameF,'.mp4" -framerate 10  -i "video_sync\\vetere\\',iNameF,'\\frame_%05d.png" -r ',frameRate,'  -pix_fmt yuva420p  -filter_complex "[1:v]format=yuva420p,colorchannelmixer=aa=0.8[ckout];[0:v][ckout]overlay[res];[res]setpts=1*PTS[res]"  -map "[res]"  -c:v libx264 -map 0:a -c:a copy "video_sync\\vetere\\_final videos\\Vetere Seduta ',iNameF,' Sync.mp4"')
  # # ffmpeg = 'cd C:/videoConversion && ffmpeg -y  -i video_sync/cc_video_deleteme/S01.mp4 -framerate 25  -i video_sync/S01/frame_%05d.png -r 25  -pix_fmt yuva420p  -filter_complex "[1:v]format=yuva420p,colorchannelmixer=aa=0.8[ckout];[0:v][ckout]overlay" -c:v libx264 video_sync/output/S01_AMICO1.1lag4.mp4'

  x = system("cmd.exe",input=ffmpeg);beepr::beep("fanfare")
 
  
  
  #
  #
  #   ffmpeg
  #   -y                                                                       | replace output file
  #   -thread_queue_size 1024                                                  | improve performance?
  #   -i "...joined.mp4"                                                       | input 0, video originale (25fps)
  #   -framerate 10  -i ".../frame_%05d.png" -r 25  -pix_fmt yuva420p          | input 1, combina [-framerate n] immagini al secondo in un video a [-r n] fps
  #   -filter_complex "[1:v]format=yuva420p,colorchannelmixer=aa=0.8[ckout];   | applica trasparenza a input1 in oggetto temp ckout
  #   [0:v][ckout]overlay[res];                                                | overlay di video originale e ckout in oggetto temp res
  #   [res]setpts=1*PTS[res]"                                                  | accelera, rallenta video. es: 0.5*PTS[res] è 2x velocità
  #   -map "[res]" -c:v libx264                                                | codec dell'output video finale
  #   -map 0:a -c:a copy                                                       | copia lo stream audio originale
  #   "... .mp4"')                                                              | file di output



  # ffmpeg
  # -i compare-emotions/seduta90.mp4  											                      	| input 0, video originale (25fps)
  # -framerate 5  -i compare-emotions/still/frame_%05d.png -r 25  -pix_fmt yuva420p | input 1, combina 5 immagini al secondo in un video a 25fps
  # -filter_complex "[1:v]format=yuva420p,colorchannelmixer=aa=0.8[ckout];	    		| applica trasparenza a input1 in oggetto temp ckout
  # 				 [0:v][ckout]overlay"					                                      		| overlay di video originale e ckout
  # -t 10																                                      			| numero di secondi di anteprima
  # -c:v libx264													                                   				| codec dell'output finale
  # output-video-test.mp4



}

