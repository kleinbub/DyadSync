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
#'
#' @param session 
#' @param signals 
#' @param inputDir 
#' @param outputDir 
#' @param inputFilePattern 
#' @param idOrder 
#' @param idSep 
#' @param overwrite logical. Should existing files be deleted? 
#' @param sync 
#' @param WIN 
#' @param w 
#' @param h 
#' @param frameRate 
#' @param layout 
#' @param col_sync 
#' @param signal_height 0 to 1. vertical space of each signal, relative to the video frame
#' @param ffmpeg.query "default" or a custom FFmpeg
#' @param ffmpeg_bin 'PATH' if FFmpeg is added to system environment. else full path to FFmpeg's bin directory.
#' @param multi.core numeric to specify a number of cores to be used, or logical to enable/disable automatic multi-core computation.
#' @param speed numeric
#' @param keepAudio logical
#' @param ... 
#' @details
#' Default FFmpeg query:
#' ffmpeg -y -thread_queue_size 1024
#' -i "inputDir/inputFile.mp4"
#' -framerate 25 -i "framesDir/framesFileFormat"
#' -r 25 -pix_fmt yuva420p
#' -filter_complex "
#'    [0:v]scale=2*trunc(iw*sar/2):ih,setsar=1[0v_1];
#'    [1:v]format=yuva420p,colorchannelmixer=aa=0.8[ckout];
#'    [0v_1][ckout]overlay[resv]; [resv]setpts=1*PTS[resv];
#'    [0:a]atempo=1[resa];"
#' -map "[resv]" -c:v libx264
#' -map "[resa]"
#' "outputDir/outputFile.mp4"

#' @return
#' @export
#' 
#' @examples

drawFrames = function(session, signals, 
                      inputDir, outputDir,
                      inputFilePattern = ".mp4",
                      idOrder=c(F,F,T), idSep="_",
                      overwrite = FALSE,
                      sync=NA, WIN=30,
                      w=1280, h=720, frameRate = "auto",
                      layout = c(1), col_sync = "default",
                      signal_height = 0.3,
                      ffmpeg.query="default", ffmpeg_bin="PATH", multi.core = TRUE, 
                      speed=1, keepAudio = TRUE, ... ){
  ############

  # session = lr[[i]];signals = "SC";
  # inputDir=inputDir; outputDir=outputDir;
  # inputFilePattern = "converted.mp4";
  # idOrder=c(T,T); idSep="_";
  # overwrite = FALSE;
  # sync="amico1"; WIN=30;
  # w=1280; h=720; frameRate = "auto";
  # layout = c(1); col_sync = "default";
  # signal_height = 0.3;
  # ffmpeg.query="default"; ffmpeg_bin="PATH"; multi.core = TRUE;
  # speed=1; keepAudio = TRUE
  #############
  
  #################################################
  ## SETUP
  
  iNameF  = UID(session)
  
  #paths THIS PROBABLY BE CHANGED
  work_dir = file.path(outputDir, iNameF) #session directory
  frames_dir = file.path(work_dir,"frames")
  #a better logic probably is: IF folder is there, skip frame generation and only do ffmpeg
  
  if(!dir.exists(work_dir)){
    dir.create(frames_dir, showWarnings = T, recursive=T)
    # stop(paste(work_dir, "was just created. Please put the video file inside it."))
  } 
  frames_exists = if(dir.exists(frames_dir) && length(list.files(frames_dir))>0) TRUE else FALSE
  if(overwrite && frames_exists){
    unlink(frames_dir, recursive = TRUE)
  }
  dir.create(frames_dir, showWarnings = F, recursive=T)
  
  #check video files
  inputVideoName = list.files(inputDir, pattern = inputFilePattern,full.names = FALSE)
  inputVideoNameParts = strsplit(inputVideoName, paste0("[\\.,",idSep,"]"))
  videoFile = c()
  for(l in seq_along(inputVideoName)){
    parts = inputVideoNameParts[[l]]
    if(
      (sessionId(session) %in% parts || lead0(sessionId(session), 3) %in% parts )
      && dyadId(session) %in% parts
    ){ 
      videoFile = c(videoFile, inputVideoName[l])
      ss = l
    }
  }
  if(length(videoFile)!=1) stop("multiple or no match for inputFilePattern in inputDir. ",UID(session)," was expected")
  
  inputVideo =  list.files(inputDir, pattern = ".mp4",full.names = TRUE)[ss]
  outputVideo = file.path(work_dir, paste0(iNameF,if(!is.na(sync)){paste0("_",sync)},".mp4"))
  
  
  #find ffmpeg
  if(ffmpeg_bin == "PATH"){
    # Read the user-specific PATH environment variable
    user_path <- Sys.getenv("PATH")
    # Split the PATH into individual directories
    path_dirs <- unlist(strsplit(user_path, ";"))
    # Look for the ffmpeg\bin folder
    ffmpeg_bin <- grep("ffmpeg", path_dirs, value = TRUE,ignore.case = TRUE)
    ffmpeg_exe = "ffmpeg"
    ffprobe_exe = file.path(ffmpeg_bin,"ffprobe.exe")
  } else {
    # check ffmpeg installation folder
    if(!file.exists(file.path(ffmpeg_bin,"ffmpeg.exe")) ||
       !file.exists(file.path(ffmpeg_bin,"ffprobe.exe"))) stop("up to date ffmpeg.exe and ffprobe.exe must be found in", ffmpeg_bin, ".")
    ffmpeg_exe = file.path(ffmpeg_bin,"ffmpeg.exe")
    ffprobe_exe = file.path(ffmpeg_bin,"ffprobe.exe")
  }
  
  #test ffmpeg
  if (!tryCatch({
    system2(ffmpeg_exe, args="-version", stdout = TRUE, stderr = TRUE)
    TRUE  # If it runs without error, return TRUE
  }, error = function(e) {
    FALSE 
  })) {
    stop("ffmpeg is not installed or not in the system's PATH.\n")
  }
  
  #test ffprobe
  if (!tryCatch({
    system2(ffprobe_exe, args="-version", stdout = TRUE, stderr = TRUE)
    TRUE  # If it runs without error, return TRUE
  }, error = function(e) {
    FALSE 
  })) {
    stop("ffprobe is not installed or not in the system's PATH.\n")
  }
  
  
  #extrapolate average video frameRate
  x = system2(ffprobe_exe, args=paste0('-show_streams "', inputVideo,'"'), stdout = TRUE)
  avg = sapply(x, grepl, pattern="avg_frame_rate" )
  avg = which(avg)
  # Now I am assuming that the first stream is always the video
  avg = x[avg[1]]
  avg = strsplit(avg, "=")[[1]][2]
  avg = strsplit(avg, "/")[[1]]
  avg = as.numeric(avg[1]) / as.numeric(avg[2])
  # if(abs(avg- frameRate)>1) stop("avg and frameRate mismatch")
  
  vidur =  sapply(x, grepl, pattern="duration=" , fixed =TRUE)
  vidur = which(vidur)
  vidur = x[vidur[1]]
  vidur = strsplit(vidur, "=")[[1]][2]
  vidur = as.numeric(vidur)
  if(is.na(vidur)) stop("Could not determine video duration :-(")
  
  if(frameRate == "auto" || frameRate == avg){
    frameRate = avg
    audio_speed_factor = 1
  } else {
    if(!is.numeric(frameRate)) stop("frameRate must be 'auto' or numeric")
    audio_speed_factor = frameRate / avg
  }
  
  
  #################################################
  ## FRAME DRAWING
  
  if(overwrite || !frames_exists){
    
    
    # Setup color scale for synchrony
    if(!is.na(sync)){
      if(length(col_sync)==1 && col_sync == "default") col_sync= c("#f70c3f","#c16107", "#878787","#878787","#878787", "#86a817", "#23c647")
      colfunc <- grDevices::colorRampPalette(colors=col_sync, bias=1)
      legendSteps =20
      colz = colfunc(legendSteps)
      bins = seq(1/legendSteps,1, by=1/legendSteps)
    }
  
    
    
    if(!all(signals %in% names(session))) stop("Specified signals not found in session ", iNameF)
    #restrict session to used signals
    s = session[signals] 
    #resample to video frameRate
    s = signalFilter(s, DyadSync::resample, newSampRate = frameRate)
    
    win = round(WIN*frameRate)
    win2 = trunc(win/2)
    
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
    nsignals = length(signals)
    sx = s  
    
    
    low_b = cumsum(rep(signal_height,nsignals)) -signal_height + 0.05
    high_b =  cumsum(rep(signal_height,nsignals)) + 0.05
    for(ii in seq_along(signals)){
      log = capture.output({
        sx[[ii]]$s1 = rescaleByWin(sx[[ii]]$s1,WIN=60*8, newMin=low_b[ii], newMax=high_b[ii], cores=F)
        sx[[ii]]$s2 = rescaleByWin(sx[[ii]]$s2,WIN=60*8, newMin=low_b[ii], newMax=high_b[ii], cores=F)
      })
      sx[[ii]]$s1 = window(sx[[ii]]$s1, start = frame_start, end = frame_end)
      sx[[ii]]$s2 = window(sx[[ii]]$s2, start = frame_start, end = frame_end)
    }
    
    if(nsignals>1) warning("sync and frames will be extracted from first signal only")
    signal = signals[1]
    if(!is.na(sync)){
      sx[[signal]][[sync]]$sync =  window(sx[[signal]][[sync]]$sync, start = frame_start, end = frame_end)
      sx[[signal]][[sync]]$lag  =  window(sx[[signal]][[sync]]$lag,  start = frame_start, end = frame_end)
    }
  
    #swinz describes the windows to be plotted
    swinz = nwin(sx[[signal]]$s1, WIN=WIN, INC = 1/frameRate, flex=FALSE, SR = frameRate,return = "all")
    showtimes = timeMaster(trunc(swinz$mid_t), out = "hour")
    swinz$show_t = showtimes
    nframes = nrow(swinz)
  
    #lag is stored in samples in amico1 so:
    if(!is.na(sync)){
      lagResolution = frequency(session[[1]]$s1) 
      sx[[signal]][[sync]]$sync = resample(sx[[signal]][[sync]]$sync, frameRate)
      sx[[signal]][[sync]]$lag  = resample(sx[[signal]][[sync]]$lag, frameRate)
      
      sy = rangeRescale(sx[[signal]][[sync]]$sync,0,1,-1,1)
      LAGs = sx[[signal]][[sync]]$lag / lagResolution
      sy[is.na(sy)] = 0
      max_lag = attr(sx[[signal]][[sync]], "lagSec") 
    }
    
  
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
      
      # for(j in 1:100){ ##ONLY FOR DEBUG, UNCOMMENT THIS TO RUN AS SIMPLE LOOP
      #draw blank frames
      
      # ## j iterates on all frames (pad + signal)
      # ## k iterates only on signals
      # if(j <= npadframes) {
      # 
      #     png( file.path(frames_dir,paste0("frame_",DyadSync::lead0(j,7),".png")),width=w, height = h, units = "px",type = "cairo")
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
      png( file.path(frames_dir,paste0("frame_",DyadSync::lead0(j,7),".png")),width=w, height = h, units = "px",type = "cairo")
      par(bg=NA,mar=c(0,0,0,0),oma=c(0,0,0,0))
      plot(-100, xlim=c(0,1),ylim=c(0,1), axes = F, frame.plot = F, xlab = "",ylab="", xaxs = "i", yaxs = "i")
      points(0.85, 0.953, cex = 3, pch = 16, col="#ff1111")
      text(0.95, 0.95, showtimes[j],cex=2,pos = 2)
      # text(0.95,0.95, pos=2, labels = "")
      
      for(ii  in 1:nsignals){
        si1 = window(sx[[ii]]$s1, start = swinz$start_t[j], end = swinz$end_t[j])
        # test = window(s$SC$s1, start = swinz$start_t[j], end = swinz$end_t[j])
        si2 = window(sx[[ii]]$s2, start = swinz$start_t[j], end = swinz$end_t[j])
        if(!is.na(sync)){
          syi = window(sy        , start = swinz$start_t[j], end = swinz$end_t[j])
          lagi= window(LAGs      , start = swinz$start_t[j], end = swinz$end_t[j])
        }
        xvals = seq(0.05,0.95,length.out = win)
        lines(xvals, si1$y, col=2,lwd=5)
        lines(xvals, si2$y, col=4,lwd=5)
        
        #il valore di lag in secondi
        if(!is.na(sync)){wlag = DyadSync::at(LAGs, time = swinz$mid_t[j])$y} else wlag = NA
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
          
        }
        
        text(0.00,0.95,paste0("~ ", DyadSync::s1Name(sx[[1]])),cex=3,col=2, pos = 4)
        text(0.00,0.90,paste0("~ ", DyadSync::s2Name(sx[[1]])),cex=3,col=4, pos = 4)
        
        abline(h = low_b)
        abline(h = high_b)
      }
      dev.off()
      
    }
    parallel::stopCluster(cl); graphics.off();beepr::beep()
    
  }
  
  complexFilterVideo = paste(
    # "[0:a]asetpts=",1/speed,"*PTS[resa];",
    "[0:v]scale=2*trunc(iw*sar/2):ih,setsar=1[0v_1];",
    # "[0:a]atempo=2.0[resa]",
    "[1:v]format=yuva420p,colorchannelmixer=aa=0.8[ckout];",
    # "[1:v]format=argb,geq=r='r(X,Y)':a='0.5*alpha(X,Y)'[ckout];",
    
    "[0v_1][ckout]overlay[resv];",
    paste0("[resv]setpts=",1/speed,"*PTS[resv];"),
    "[0:a]atempo=",speed,",atempo=",audio_speed_factor,"[resa];",#"[resa]pan=stereo|FL=FL|FR=FR[resa];",
    "")
  
  
  
  # if(keepAudio){
  #   ffmpeg.audio = "-map 0:a -c:a copy" 
  # } else {
  #   ffmpeg.audio = ""
  # }
  # 
  
  
  ffmpeg.query = paste0(
    if(overwrite){'-y '},
    '-thread_queue_size 1024 ',
    '-i "',inputVideo, '" ',
    '-framerate ',frameRate,' ', '-i "',frames_dir,'\\frame_%07d.png" ',
    '-r ',frameRate,' ','-pix_fmt yuva420p ',
    # '-f lavfi -t 0.1 -i anullsrc',
    '-filter_complex "',complexFilterVideo,'"  ',
    # '-filter_complex "',complexFilterAudio,'"  ',
    '-map "[resv]" -c:v libx264 ',
    '-map "[resa]" ',
    '"',outputVideo,'"')
  
  x = system2(ffmpeg_exe, args = ffmpeg.query)
  # x = system("cmd.exe",input=ffmpeg);
  beepr::beep("fanfare")

  
}

