#!/usr/bin/env Rscript

# Tracking and characterisation of pre-identified particles

posttrack = function(remove_drift = FALSE){
  
  # source all scripts in the current directory
  setwd("C:\\Users\\Peter\\Documents\\Uni\\PhD\\2d-particle-tracking")
  filelist <- c("characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "characterisation/msd.r", "characterisation/shift.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  varfilename <- "C:\\Users\\Peter\\Dropbox\\tracking\\averaged_coords.txt"
  #put slash on end of dirname
  vardirname <- "C:\\Users\\Peter\\Dropbox\\tracking\\"
  
  #Track variables
  varedgecutoff <- 10     #Set to same value as used in particle identification
  varmaxdisp <- 5         #Used in tracking - the maximum allowed interframe displacement
  
  #Other variables that I can't think of a title for
  varparticlesize = 15    #Used as the wavevector for isf
  vartimestep = 1     #Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  vargofrframes = 10      #How many frames of data to analyse for the g(r)
  
  
  ### Main ###
  
  # Particles come in three tab seperated columns: x, y, framenumber
  pt <- read.delim(varfilename, header=FALSE, sep="\t")
  # Bulk out pt with a few blank columns so that iantrack doesn't complain
  pt <- cbind(pt[, 1:2], matrix(0,nrow = length(pt[, 1]), ncol = 3), pt[, 3])
  
  varimages <- max(pt[, 6])
  varimgx <- max(pt[, 1])          #Width of the image in pixels
  varimgy <- max(pt[, 2])          #Height of the image in pixels
  
  
  #Get particle count for each image
  particlecount <- rep(0, varimages-1)
  for(i in 0:varimages-1) {particlecount[i] <- sum(pt[, 6] == i)}
  particlecount <- cbind(seq(1,varimages-1), seq(1,varimages-1)*vartimestep, matrix(particlecount))
  particlecount <- rbind(c("Frames", "Time", "Particle Count"), particlecount)
  write(t(particlecount),file=paste(vardirname, "particlecount.txt", sep=""),ncolumns=3,sep="\t")
  
  tr <- iantrack(pretrack=pt[],maxdisp=varmaxdisp,imgsize=c(varimgx+varedgecutoff,varimgy+varedgecutoff),goodenough=3)
  if(remove_drift == TRUE) 
  {trnodrift <- driftremoval(tr, display_graphs = FALSE)}
  
  write(t(tr),file=paste(vardirname, "tracked_coords.txt", sep=""),ncolumns=7,sep="\t")
  #Don't need to filter tr as it was already filtered from pt
  #Do msd and write it out
  
  fulltrajectoriesdata = matrix(data = 0, nrow = 1, ncol = 7)
  
  # Get indicies of particles in first frame
  browser()
  initialparticles = tr[which(tr[, 6] == 0), 7]
  finalparticles = tr[which(tr[, 6] == varimages), 7]
  fulltrajectories = intersect(initialparticles, finalparticles)
  for(i in 1:length(fulltrajectories))
    fulltrajectoriesdata = rbind(fulltrajectoriesdata, tr[which(tr[, 7] == fulltrajectories[i]),])
  
  #remove blank first row
  fulltrajectoriesdata = fulltrajectoriesdata[-1,]
    
  if(remove_drift == TRUE) 
  {msq <- msd(trnodrift)}
  else
  {msq <- msd(tr)}
  msq <- cbind(msq[,1], msq[,1]*vartimestep, msq[,2:6], msq[,6]/varparticlesize, msq[,7]) # multiply the first column by vartimestep to give time in seconds, divide the 6th column to give msd in dimeters
  msq <- rbind(c("Frames", "Time", "Mean dx", "Mean dy", "Mean dx^2", "Mean dy^2", "Mean dr^2 (Pixels)", "Mean dr^2 (Diameters)", "Particle Count"), msq)
  write(t(msq),file=paste(vardirname, "msd.txt", sep=""),ncolumns=9,sep="\t")
  
  #Do isf and write it out
  if(remove_drift == TRUE) 
  {fsqt <- isf(trnodrift,length=varparticlesize)}
  else
  {fsqt <- isf(tr,length=varparticlesize)}
  fsqt <- cbind(fsqt[,1], fsqt[,1]*vartimestep, fsqt[,2:5]) # multiply the first column by vartimestep to give time in seconds
  fsqt <- rbind(c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"), fsqt)
  write(t(fsqt), file=paste(vardirname, "isf.txt", sep=""), ncolumns=6, sep="\t")
  
  #Do g(r) and write it out
  gr <- gr2d(tr,nbins=200,deltar=0.5,imgsize=c(varedgecutoff,varedgecutoff,varimgx,varimgy), nframes=vargofrframes)
  gr <- cbind(gr[,1], gr[,1]/varparticlesize, gr[,2]) # multiply the first column by vartimestep to give time in seconds
  gr <- rbind(c("Distance (Diameters)", "Distance (Pixels)", "g(r)"), gr)
  write(t(gr),file=paste(vardirname, "gr.txt", sep=""),ncolumns=3,sep="\t")
  
}
