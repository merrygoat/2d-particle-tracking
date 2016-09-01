#!/usr/bin/env Rscript

# Characterisation of pre-identified and pre-tracked particles

postposttrack = function(){
  
  # source all scripts in the current directory
  #setwd("C:\\Users\\Peter\\Documents\\Uni\\PhD\\2d-particle-tracking")
  setwd("/Users/pc9836/Documents/git/2d-particle-tracking")
  filelist <- c("characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "characterisation/msd.r", "characterisation/shift.r")
  sapply(filelist,source,.GlobalEnv)
  
  #File directory variables
  #varfilename <- "C:\\Users\\Peter\\Dropbox\\tracking\\averposaged_coords.txt"
  varfilename <- "/Users/pc9836/Documents/LAMMPS/simulations/hard_spheres/55/6/6.0depthslice.xmol"
  #put slash on end of dirname
  #vardirname <- "C:\\Users\\Peter\\Dropbox\\tracking\\"
  vardirname <- "/Users/pc9836/Documents/LAMMPS/simulations/hard_spheres/55/6/"

  #Track variables
  varedgecutoff <- 0     #Set to same value as used in particle identification
  varmaxdisp <- 0.3         #Used in tracking - the maximum allowed interframe displacement
  
  #Other variables that I can't think of a title for
  varparticlesize = 1  #Used as the wavevector for isf
  vartimestep = 1     #Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  vargofrframes = 5      #How many frames of data to analyse for the g(r)
  
  
  ### Main ###
  
  # Particles come in five tab seperated columns: x, y, z, framenumber, particlenumber
  pt <- read.delim(varfilename, header=FALSE, sep=" ")
  # Bulk out pt with a few blank columns so that iantrack doesn't complain
  pt <- cbind(pt[, 1:2], matrix(0,nrow = length(pt[, 1]), ncol = 3), pt[, 4:5])
  
  varimages <- max(pt[, 6])
  varimgx <- max(pt[, 1])          #Width of the image in pixels
  varimgy <- max(pt[, 2])          #Height of the image in pixels
  
  #Get particle count for each image
  particlecount <- rep(0, varimages-1)
  for(i in 0:varimages-1) {particlecount[i] <- sum(pt[, 6] == i)}
  particlecount <- cbind(seq(1,varimages-1), seq(1,varimages-1)*vartimestep, matrix(particlecount))
  particlecount <- rbind(c("Frames", "Time", "Particle Count"), particlecount)
  write(t(particlecount),file=paste(vardirname, "particlecount.txt", sep=""),ncolumns=3,sep="\t")
  
  fulltrajectoriesdata = matrix(data = 0, nrow = 1, ncol = 7)
  
  #msq <- msd(pt)
  #msq <- cbind(msq[,1], msq[,1]*vartimestep, msq[,2:6], msq[,6]/varparticlesize, msq[,7]) # multiply the first column by vartimestep to give time in seconds, divide the 6th column to give msd in dimeters
  #msq <- rbind(c("Frames", "Time", "Mean dx", "Mean dy", "Mean dx^2", "Mean dy^2", "Mean dr^2 (Pixels)", "Mean dr^2 (Diameters)", "Particle Count"), msq)
  #write(t(msq),file=paste(vardirname, "msd.txt", sep=""),ncolumns=9,sep="\t")
  
  fsqt <- isf(pt,length=varparticlesize)
  fsqt <- cbind(fsqt[,1], fsqt[,1]*vartimestep, fsqt[,2:5]) # multiply the first column by vartimestep to give time in seconds
  fsqt <- rbind(c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"), fsqt)
  write(t(fsqt), file=paste(vardirname, "isf.txt", sep=""), ncolumns=6, sep="\t")
  
  #Do g(r) and write it out
  #gr <- gr2d(tr,nbins=200,deltar=0.05,imgsize=c(varedgecutoff,varedgecutoff,varimgx,varimgy), nframes=vargofrframes)
  #gr <- cbind(gr[,1], gr[,1]/varparticlesize, gr[,2]) # multiply the first column by vartimestep to give time in seconds
  #gr <- rbind(c("Distance (Diameters)", "Distance (Pixels)", "g(r)"), gr)
  #write(t(gr),file=paste(vardirname, "gr.txt", sep=""),ncolumns=3,sep="\t")
  
}
