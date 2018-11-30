#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

confocaltrackroutine = function(remove_drift = TRUE){
  
  base_directory <- toString(strsplit(getSrcDirectory(function(x) {x}), "/confocal"))
  cat(base_directory)
  # source all scripts in the current directory
  setwd(base_directory)
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r", "characterisation/psi6series.r", "characterisation/psi6loc.r", "characterisation/drawvoronoi.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)

  #File directory variables
  varfilename <- "./documentation/t_"
  
  #put slash on end of dirname
  vardirname <- "./sample/"
  
  #Pretrack variables
  varimages <- 821        #How many image to read from varfilename
  vardiameter <- 11       #Particle diameter - used in particle identification
  varfilter <- 11         #Parameter for lowpass filter
  varbgavg <- 11           #Parameter for lowpass filter
  varmasscut <- 3         #Lowest integrated brightness for particle
  varminimum <- 0.1       #Lowest pixel value for center of particle
  
  #Track variables
  varedgecutoff <- 20     #Cuts off this many pixels from each edge of the image in all data output - this is because particle identification is bad around the edges.
  varmaxdisp <- 800         #Used in tracking - the maximum allowed interframe displacement
  
  #Other variables that I can't think of a title for
  varparticlesize = 11    #Used as the wavevector for isf
  vartimestep = 0.0294 #0.05656    #Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  vargofrframes = 821     #How many frames of data to analyse for the g(r)
  
  
  ### Main ###
  
  testimg <- readImage(paste(varfilename,"00000.tif",sep=""))
  
  varimgx <- dim(testimg)[1]          #Width of the image in pixels
  varimgy <- dim(testimg)[2]          #Height of the image in pixels
  
  # Identify particles in images
  pt <- pretrack(filename=varfilename,images=varimages,diameter=vardiameter,filter=varfilter,bgavg=varbgavg, masscut=varmasscut,minimum=varminimum,chan="grey")
  ptfilt <- which(pt[,1] > varedgecutoff & pt[,1] < (varimgx-varedgecutoff) & pt[,2] > varedgecutoff & pt[,2] < (varimgy-varedgecutoff))
  
  #Print out some sample overcircled images so the user can check the tracking is OK.  
  
  #Remove /10*varimages to write all the overcicled images. Change 9 for vargofrframes in loop
  for(i in 0:9) {
    imagecn <- paste(varfilename,formatC(round(i/10*varimages),flag="0",digits=4),".tif",sep="") 
    imagecn <- channel(readImage(imagecn), "gray")
  
    #r <- channel(readImage(imagecn), "r")
    #g <- channel(readImage(imagecn), "g")
    #b <- channel(readImage(imagecn), "b")
    
    #imagecn <- 0.21*r+0.71*g+0.07*b
    
    temp <- pt[ptfilt,]
    ptmask <- which(temp[,6] == round(i/10*varimages))
    circledimage <- overcirc(imagecn, temp[ptmask,], rad=5) 
    #writeImage(circledimage, paste(vardirname, "overcirc", formatC(round(i/10*varimages),flag="0",digits=4), ".tif",sep=""))
  }
  
  
  #Get particle count for each image
  particlecount <- rep(0, varimages-1)
  for(i in 0:varimages-1) {particlecount[i] <- sum(pt[ptfilt,6] == i)}
  particlecount <- cbind(seq(1,varimages-1), seq(1,varimages-1)*vartimestep, matrix(particlecount))
  write.table(particlecount, file=paste(vardirname, "particlecount.txt", sep=""), row.names=FALSE, col.names=c("Frames", "Time", "Particle Count"))
  
  # Track trajectories of identified particles
  tr <- iantrack(pretrack=pt[ptfilt,], maxdisp=varmaxdisp, imgsize=c(varimgx,varimgy), goodenough=5)
  
  if(remove_drift == TRUE) {
    trnodrift <- driftremoval(tr, display_graphs = FALSE)
  }
  
  write.table(tr, file=paste(vardirname, "raw_coords.txt", sep=""), row.names=FALSE, col.names=FALSE)
  
  #Don't need to filter tr as it was already filtered from pt
  
  #Do msd and write it out
  if(remove_drift == TRUE)
  {msq <- msd(trnodrift)}
  else
  {msq <- msd(tr)}
  msq <- cbind(msq[,1], msq[,1]*vartimestep, msq[,2:6], msq[,6]/varparticlesize, msq[,7]) # multiply the first column by vartimestep to give time in seconds, divide the 10th column to give msd in dimeters
  write.table(msq, file=paste(vardirname, "msd.txt", sep=""), row.names=FALSE, 
              col.names=c("Frames", "Time", "Mean dx", "Mean dy", "Mean dx^2", "Mean dy^2", "Mean dr^2 (Pixels)", "Mean dr^2 (Diameters)", "Particle Count"))

  
  #Do isf and write it out
  if(remove_drift == TRUE) 
  {fsqt <- isf(trnodrift,length=19)}
  else
  {fsqt <- isf(tr,length=19)}
  fsqt <- cbind(fsqt[,1], fsqt[,1]*vartimestep, fsqt[,2:5]) # multiply the first column by vartimestep to give time in seconds
  
  write.table(fsqt, file=paste(vardirname, "isf.txt", sep=""), row.names=FALSE, col.names=c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"))
  

  #Do g(r) and write it out
  gr <- gr2d(tr,nbins=200, deltar=0.5, imgsize=c(varedgecutoff, varedgecutoff, varimgx-varedgecutoff, varimgy-varedgecutoff), nframes=vargofrframes)
  gr <- cbind(gr[,1], gr[,1]/varparticlesize, gr[,2]) # multiply the first column by vartimestep to give time in seconds
  write.table(gr, file=paste(vardirname, "gr.txt", sep=""), row.names=FALSE, col.names=c("Distance (Pixels)", "Distance (Diamaters)", "g(r)"))
  
  #Do psi_6
  psi6 <- psi6series(tr)
  write.table(psi6, file=paste(vardirname, "psi6.txt", sep=""), row.names=FALSE, 
              col.names=c("xcoords", "ycoords", "psi_6 real", "psi_6 imaginary", "psi_6 modulus", "frame number"))
  
  #Print the modulus
  cat ("The mean modulus of psi6 is: ", mean(psi6[2:nrow(psi6),5]), "\n")

  #Voronoi
  #vor <- drawvoronoi(psi6[2:nrow(psi6), 1], psi6[2:nrow(psi6), 2], psi6[2:nrow(psi6), 5], imgsize=c(0,varimgx,0,varimgy))
  #imgagev <- vor 
  #writeImage(imagev,file=paste(vardirname,"voronoi"), ".tif", sep="")
  
}
