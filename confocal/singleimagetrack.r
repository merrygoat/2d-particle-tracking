#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

confocaltrackroutine = function(remove_drift = TRUE){
  
  # source all scripts in the current directory
  setwd("/Users/mm8769/Documents/Analysis/RParticletracking")
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  varfilename <- "/Users/mm8769/Documents/Images/Microfluidics/2016_04_12/open_init_8_2.tif"
  #put slash on end of dirname
  vardirname <- "/Users/mm8769/Documents/Images/Microfluidics/2016_04_12/"
  
  ### Main ###
  
  img <- channel(readImage(varfilename), "grey")
  
  varimgx <- dim(img)[1]          #Width of the image in pixels
  varimgy <- dim(img)[2]          #Height of the image in pixels
  
  lpbf <- lowpass(img,lobject=19,bgavg=19)
  #display(lpbf)

  fbf <- feature(lpbf,diameter=20,masscut=3,minimum=0.3,verbose = TRUE)
  
  display(overcirc(lpbf,fbf,rad=10))
  #browser()
  tmp <- matrix(data = 0, ncol = 1, nrow = dim(fbf)[1])
  
  fbf <- cbind(fbf,tmp)            
  
  #Do g(r) and write it out
  gr <- gr2d(fbf,nbins=150,deltar=0.5,imgsize=c(0,0,varimgx,varimgy), nframes=1)
  gr <- rbind(c("Distance (Pixels)", "g(r)"), gr)
  write(t(gr),file=paste(vardirname, "gr.txt", sep=""),ncolumns=2,sep="\t")

}
