#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

singletrackroutine = function(){
  
  base_directory <- toString(strsplit(getSrcDirectory(function(x) {x}), "/confocal"))
  cat(base_directory)
  # source all scripts in the current directory
  setwd(base_directory)
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r", "characterisation/psi6series.r", "characterisation/psi6loc.r", "characterisation/drawvoronoi.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)

  #File directory variables
  varfilename <- "./documentation/brightfield.tif"
  #put slash on end of dirname
  vardirname <- "./sample/"
  
  vardiameter <- 9                # Particle size in pixels
  
  ### Main ###
  
  img <- channel(readImage(varfilename), "grey")
  varimgx <- dim(img)[1]          # Width of the image in pixels
  varimgy <- dim(img)[2]          # Height of the image in pixels
  
  # Apply the lowpass filter to blur the image
  lpbf <- lowpass(img,lobject=11,bgavg=11)

  # Identify particles from the blurred image
  fbf <- feature(lpbf, diameter=9, masscut=5, minimum=0.1, verbose = TRUE)
  
  
  display(overcirc(lpbf, fbf, rad=vardiameter/2))
  
  # Add a blank column to the tracked particle data so it the correct shape for the analysis functions
  fbf <- cbind(fbf, matrix(data = 0, ncol = 1, nrow = dim(fbf)[1])) 
  
  #Do g(r) and write it to a file
  gr <- gr2d(fbf, nbins=150, deltar=0.5, imgsize=c(0,0,varimgx,varimgy), nframes=1)
  gr <- cbind(gr[,1], gr[,1]/vardiameter, gr[,2])
  write.table(gr, file=paste(vardirname, "gr.txt", sep=""), row.names=FALSE,
              col.names=c("Distance (Pixels)", "Diastance (Diameters)",  "g(r)"))
  
  #Do psi_6
  psi6 <- psi6loc(fbf[,1],fbf[,2])
  write.table(psi6, file=paste(vardirname, "psi6.txt", sep=""), row.names=FALSE,
              col.names=c("xcoord", "ycoord", "real psi_6", "imaginary psi_6", "modulus"))
  
  #Print the mean value of psi_6:
  cat ("The mean modulus of psi6 is:", mean(psi6[, 5]), "\n")

  #Voronoi
  #vor <- drawvoronoi(psi6[2:nrow(psi6), 1], psi6[2:nrow(psi6), 2], psi6[2:nrow(psi6), 5], imgsize=c(0,varimgx,0,varimgy))

}