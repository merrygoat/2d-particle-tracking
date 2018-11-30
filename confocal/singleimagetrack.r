#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

singletrackroutine = function(remove_drift = TRUE){
  
  # source all scripts in the current directory
  setwd("/Users/am15146/Documents/Tracking_codes/2d-particle-tracking/")
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r", "characterisation/psi6series.r", "characterisation/psi6loc.r", "characterisation/drawvoronoi.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  library(tripack)
  
  #File directory variables
  varfilename <- "/Users/am15146/Documents/Research/PMMA/2016/September_2016/13_09_16/Mov_2/clusters/t_00000.tif"
  #put slash on end of dirname
  vardirname <- "/Users/am15146/Documents/Research/PMMA/2016/September_2016/13_09_16/Mov_2/clusters/"
  
  ### Main ###
  
  img <- channel(readImage(varfilename), "grey")
                 
  vardiameter <- 11               #Particle size in pixels
  
  varimgx <- dim(img)[1]          #Width of the image in pixels
  varimgy <- dim(img)[2]          #Height of the image in pixels
  
  
  lpbf <- lowpass(img,lobject=11,bgavg=11)
  

  #display(lpbf)

  fbf <- feature(lpbf,diameter=11,masscut=5,minimum=0.1,verbose = TRUE)
  
  
  display(overcirc(lpbf,fbf,rad=6))
  tmp <- matrix(data = 0, ncol = 1, nrow = dim(fbf)[1])
  
  fbf <- cbind(fbf,tmp) 
  
  #Do g(r) and write it out
  gr <- gr2d(fbf,nbins=150,deltar=0.5,imgsize=c(0,0,varimgx,varimgy), nframes=1)
  gr <- cbind(gr[,1], gr[,1]/vardiameter, gr[,2]) # multiply the first column by vartimestep to give time in seconds
  gr <- rbind(c("Distance (Pixels)", "Diastance (Diameters)",  "g(r)"), gr)
  #write(t(gr),file=paste(vardirname, "gr.txt", sep=""),ncolumns=3, sep="\t")
  
  #Do psi_6
  psi6 <- psi6loc(fbf[,1],fbf[,2])
  psi6 <- rbind(c("xcoord","ycoord", "psi_6 real", "psi_6 imaginary", "psi_6 modulus"), psi6)
  write(t(psi6),file=paste(vardirname, "psi6.txt", sep=""),ncolumns=5,sep="\t")
  #Compute the mean value of psi_6:
  # 1- convert to numbers!
  class(psi6)<-"numeric"
  #Take the modulus
  cat ("The mean modulus of psi6 is:")
  cat (mean(psi6[2:nrow(psi6),5]))
  #return(psi6)
  
  #Voronoi
  #vor <- drawvoronoi(psi6[2:nrow(psi6), 1], psi6[2:nrow(psi6), 2], psi6[2:nrow(psi6), 5], imgsize=c(0,varimgx,0,varimgy))
  

  
}