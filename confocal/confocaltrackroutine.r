#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

confocaltrackroutine = function(remove_drift = TRUE){
  
  # source all scripts in the current directory
  setwd("/Users/am15146/Documents/Tracking_codes/2d-particle-tracking/")
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r", "characterisation/psi6series.r", "characterisation/psi6loc.r", "characterisation/drawvoronoi.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  varfilename <- "/Users/am15146/Documents/Research/PMMA/2016/December_2016/10_12_16/16.2um_gapcell/0.1/Data_analysis/13_crystal_amoeba/t_"
  
  #put slash on end of dirname
  vardirname <- "/Users/am15146/Documents/Research/PMMA/2016/December_2016/10_12_16/16.2um_gapcell/0.1/Data_analysis/13_crystal_amoeba/data/"
  
  #Pretrack variables
  varimages <- 821    #How many image to read from varfilename
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
  particlecount <- rbind(c("Frames", "Time", "Particle Count"), particlecount)
  #write(t(particlecount),file=paste(vardirname, "particlecount.txt", sep=""),ncolumns=3,sep="\t")
  
  tr <- iantrack(pretrack=pt[ptfilt,],maxdisp=varmaxdisp,imgsize=c(varimgx,varimgy),goodenough=5)
  
  
  if(remove_drift == TRUE) 
  {trnodrift <- driftremoval(tr, display_graphs = FALSE)}
  
  
  write(t(tr),file=paste(vardirname, "raw_coords.txt", sep=""),ncolumns=7,sep=" ")
  #Don't need to filter tr as it was already filtered from pt
  
  #Do msd and write it out
  if(remove_drift == TRUE)
  {msq <- msd(trnodrift)}
  else
  {msq <- msd(tr)}
  msq <- cbind(msq[,1], msq[,1]*vartimestep, msq[,2:6], msq[,6]/varparticlesize, msq[,7]) # multiply the first column by vartimestep to give time in seconds, divide the 10th column to give msd in dimeters
  msq <- rbind(c("Frames", "Time", "Mean dx", "Mean dy", "Mean dx^2", "Mean dy^2", "Mean dr^2 (Pixels)", "Mean dr^2 (Diameters)", "Particle Count"), msq)
  write(t(msq),file=paste(vardirname, "msd_.txt", sep=""),ncolumns=9,sep="\t")
  
  
  #Do isf and write it out
  if(remove_drift == TRUE) 
  {fsqt <- isf(trnodrift,length=19)}
  else
  {fsqt <- isf(tr,length=19)}
  fsqt <- cbind(fsqt[,1], fsqt[,1]*vartimestep, fsqt[,2:5]) # multiply the first column by vartimestep to give time in seconds
  fsqt <- rbind(c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"), fsqt)
  write(t(fsqt), file=paste(vardirname, "isf.txt", sep=""), ncolumns=6, sep="\t")

  
  #Do g(r) and write it out
  gr <- gr2d(tr,nbins=200,deltar=0.5,imgsize=c(varedgecutoff,varedgecutoff,varimgx-varedgecutoff,varimgy-varedgecutoff), nframes=vargofrframes)
  gr <- cbind(gr[,1], gr[,1]/varparticlesize, gr[,2]) # multiply the first column by vartimestep to give time in seconds
  gr <- rbind(c("Distance (Pixels)", "Distance (Diamaters)", "g(r)"), gr)
  write(t(gr),file=paste(vardirname, "gr.txt", sep=""),ncolumns=3,sep="\t")
  
  #Do psi_6
  psi6 <- psi6series(tr)
  psi6 <- rbind(c("xcoords", "ycoords", "psi_6 real", "psi_6 imaginary", "psi_6 modulus", "frame number"), psi6)
  write(t(psi6),file=paste(vardirname, "psi6.txt", sep=""),ncolumns=6,sep="\t")
  #Compute the mean value of psi_6:
  # 1- convert to numbers!
  class(psi6)<-"numeric"
  #Take the modulus
  cat ("The mean modulus of psi6 is: ")
  cat (mean(psi6[2:nrow(psi6),5]))
  
  #Voronoi
  #vor <- drawvoronoi(psi6[2:nrow(psi6), 1], psi6[2:nrow(psi6), 2], psi6[2:nrow(psi6), 5], imgsize=c(0,varimgx,0,varimgy))
  #imgagev <- vor 
  #writeImage(imagev,file=paste(vardirname,"voronoi"), ".tif", sep="")
  
  #Measure the quality of our tracking by dividing dyanmic samples by msd
  trackingquality <- matrix(ncol=5,nrow=varimages-1)
  trackingquality[,1:2] <- msq[-1,8:9]
  trackingquality[,3] <- as.numeric(trackingquality[,2])/as.numeric(trackingquality[1,2])
  trackingquality[,4] <- abs(as.numeric(trackingquality[,3])-(1/exp(1)))
  cat("       Diameters diffusion by 1/e complete particle tracks = ", round(as.numeric(trackingquality[which(trackingquality[,4] == min(trackingquality[,4])),1]),2))
  
}
