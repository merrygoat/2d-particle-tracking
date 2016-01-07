#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

threeDconfocaltrackroutine = function(remove_drift = TRUE){
  
  # source all scripts in the current directory
  setwd("/Users/pc9836/Documents/git/2d-particle-tracking")
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  varfilename <- "/Volumes/WIN_DATA/Confocal/STED/15-12-21/images/FITC 19"
  #put slash on end of dirname
  vardirname <- "/Volumes/WIN_DATA/Confocal/STED/15-12-21/"
  
  #Pretrack variables
  varimages <- 30        #How many image to read from varfilename
  varzdepth <- 18
  vardiameter <- 11       #Particle diameter - used in particle identification
  varfilter <- 11         #Parameter for lowpass filter
  varbgavg <- 11          #Parameter for lowpass filter
  varmasscut <- 1         #Lowest integrated brightness for particle
  varminimum <- 0.1       #Lowest pixel value for center of particle
  
  #Track variables
  varedgecutoff <- 10     #Cuts off this many pixels from each edge of the image in all data output - this is because particle identification is bad around the edges.
  varmaxdisp <- 5         #Used in tracking - the maximum allowed interframe displacement
  
  #Other variables that I can't think of a title for
  varparticlesize = 19    #Used as the wavevector for isf
  vartimestep = 60    #Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  vargofrframes = 64      #How many frames of data to analyse for the g(r)
  
  
  ### Main ###
  
  testimg <- readImage(paste(varfilename,"_t0001_z0001.png",sep=""))
  
  varimgx <- dim(testimg)[1]          #Width of the image in pixels
  varimgy <- dim(testimg)[2]          #Height of the image in pixels
  
  particlecount <- rep(0, varimages-1)
  allfsqt <- seq(from=1, to=varimages-1)
  allfsqt <- cbind(matrix(allfsqt), allfsqt*vartimestep)
  
  for(i in 1:varzdepth){
    pt <- pretrack(filename=varfilename,images=varimages,diameter=vardiameter,filter=varfilter,bgavg=varbgavg, masscut=varmasscut,minimum=varminimum,chan="grey",zstack=TRUE,znumber=i)
    ptfilt <- which(pt[,1] > varedgecutoff & pt[,1] < (varimgx-varedgecutoff) & pt[,2] > varedgecutoff & pt[,2] < (varimgy-varedgecutoff))
  
    
  #Get particle count for each image
  newparticlecount <- rep(0, varimages-1)
  for(i in 0:varimages-1) {newparticlecount[i] <- sum(pt[ptfilt,6] == i)}
  particlecount <- cbind(matrix(particlecount), matrix(newparticlecount))
    
  tr <- iantrack(pretrack=pt[ptfilt,],maxdisp=varmaxdisp,imgsize=c(varimgx,varimgy),goodenough=10)
  if(remove_drift == TRUE) 
    {trnodrift <- driftremoval(tr, display_graphs = TRUE)}
  
  #write(t(tr),file=paste(vardirname, "raw_coords.txt", sep=""),ncolumns=7,sep="\t")
  #Don't need to filter tr as it was already filtered from pt
  #Do msd and write it out
  #if(remove_drift == TRUE) 
  #  {msq <- msd(trnodrift)}
  #else
  #  {msq <- msd(tr)}
  #msq <- cbind(msq[,1], msq[,1]*vartimestep, msq[,2:6], msq[,6]/varparticlesize, msq[,7]) # multiply the first column by vartimestep to give time in seconds, divide the 6th column to give msd in dimeters
  #msq <- rbind(c("Frames", "Time", "Mean dx", "Mean dy", "Mean dx^2", "Mean dy^2", "Mean dr^2 (Pixels)", "Mean dr^2 (Diameters)", "Particle Count"), msq)
  #write(t(msq),file=paste(vardirname, "msd.txt", sep=""),ncolumns=9,sep="\t")
  
  
  #Do isf and write it out
  if(remove_drift == TRUE) 
    {fsqt <- isf(trnodrift,length=19)}
  else
    {fsqt <- isf(tr,length=19)}
  allfsqt <- cbind(matrix(allfsqt), fsqt[,4]) # multiply the first column by vartimestep to give time in seconds
  #fsqt <- rbind(c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"), fsqt)
  }
  
  write(t(fsqt), file=paste(vardirname, "isf.txt", sep=""), ncolumns=6, sep="\t")
  
  particlecount <- cbind(seq(1,varimages-1), seq(1,varimages-1)*vartimestep, matrix(particlecount))
  particlecount <- rbind(c("Frames", "Time", "Particle Count"), particlecount)
  write(t(particlecount),file=paste(vardirname, "particlecount.txt", sep=""),ncolumns=varzdepth+2,sep="\t")
  
}
