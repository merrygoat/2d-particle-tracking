#!/usr/bin/env Rscript

trackroutine = function(){
  
  # source all scripts in the current directory
  setwd("/Users/pc9836/Documents/git/2d-particle-tracking")
  filelist <- c("chi4.r", "feature.r", "gr2d.r", "iantrack.r", "isf.r", "lowpass.r", "msd.r", "pretrack.r", "shift.r", "overcirc.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  istiffstack <- TRUE     #Set to true if the images are read in as a single compound tiff image
  varfilename <- "/Volumes/WIN_DATA/Confocal/STED/15-08-04/pos10/FITC-0.52-pos10"
  vardirname <- "/Volumes/WIN_DATA/Confocal/STED/15-08-04/pos10/"
  
  #Pretrack variables
  varimages <- 200        #How many image to read from varfilename
  vardiameter <-19        #Particle diamter - used in particle identification
  varfilter <-11          #Parameter for lowpass filter
  varbgavg <- 11          #Parameter for lowpass filter
  varmasscut <- 1         #Lowest integrated brightness for particle
  varminimum <- 0.1       #Lowest pixel value for center of particle
  
  #Track variables
  varedgecutoff <- 15     #Cuts off this many pixels from each edge of the image in all data output - this is because particle identification is bad around the edges.
  varmaxdisp <- 8         #Used in tracking - the maximum allowed interframe displacement
  
  #Other variables that I can't think of a title for
  varparticlesize = 19    #Used as the wavevector for isf
  vartimestep = 1     #Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  vargofrframes = 20      #How many frames of data to analyse for the g(r)
  
  ### Main ###
  
  if (istiffstack == TRUE) {
    cat("Reading image file.\n")
    imgtiffstack <- suppressWarnings(readImage(paste(varfilename,".tif",sep="")))   #Supress to avoid warnings about unreadable metadata
  }
  else {
    imgtiffstack <- paste(varfilename,"000.tif",sep="") #If no tiff stack, read in one image anyway to get the dimensions.
  }
  
  varimgx <- dim(imgtiffstack[1])          #Width of the image in pixels
  varimgy <- dim(imgtiffstack[2])          #Height of the image in pixels
  
  pt <- pretrack(filename=varfilename,images=varimages,diameter=vardiameter,filter=varfilter,bgavg=varbgavg, masscut=varmasscut,minimum=varminimum,chan="grey", istiffstack = istiffstack, imgtiffstack = imgtiffstack)
  ptfilt <- which(pt[,1] > varedgecutoff & pt[,1] < (varimgx-varedgecutoff) & pt[,2] > varedgecutoff & pt[,2] < (varimgy-varedgecutoff))
  
  #Print out some sample overcircled images so the user can check the tracking is OK.  
  for(i in 0:9) {
    if (istiffstack == FALSE) {
      imagecn <- paste(varfilename,formatC(round(i/10*varimages),flag="0",digits=3),".tif",sep="") 
      imagecn <- channel(readImage(imagecn), "grey")
    }
    else {
      imagecn <- getFrame(imgtiffstack, round(i/10*varimages)+1)
    }
    
    temp <- pt[ptfilt,]
    ptmask <- which(temp[,6] == round(i/10*varimages))
    circledimage <- overcirc(imagecn, temp[ptmask,], rad=5) 
    writeImage(circledimage, paste(vardirname, "overcirc", formatC(round(i/10*varimages),flag="0",digits=3), ".tif",sep=""))
  }
  
  
  #Get particle count for each image
  particlecount <- rep(0, varimages-1)
  for(i in 0:varimages-1) {particlecount[i] <- sum(pt[ptfilt,6] == i)}
  particlecount <- cbind(seq(1,varimages-1), seq(1,varimages-1)*vartimestep, matrix(particlecount))
  particlecount <- rbind(c("Frames", "Time", "Particle Count"), particlecount)
  write(t(particlecount),file=paste(vardirname, "particlecount.txt"),ncolumns=3,sep="\t")
  
  tr <- iantrack(pretrack=pt[ptfilt,],maxdisp=varmaxdisp,imgsize=c(varimgx,varimgy),goodenough=10)
  
  #Don't need to filter tr as it was already filtered from pt
  #Do msd and write it out
  msq <- msd(tr)
  msq <- cbind(msq[,1], msq[,1]*vartimestep, msq[,2:6], msq[,6]/varparticlesize, msq[,7]) # multiply the first column by vartimestep to give time in seconds, divide the 6th column to give msd in dimeters
  msq <- rbind(c("Frames", "Time", "Mean dx", "Mean dy", "Mean dx^2", "Mean dy^2", "Mean dr^2 (Pixels)", "Mean dr^2 (Diameters)", "Particle Count"), msq)
  write(t(msq),file=paste(vardirname, "msd.txt"),ncolumns=9,sep="\t")
  
  #Do isf and write it out
  fsqt <- isf(tr,length=19)
  fsqt <- cbind(fsqt[,1], fsqt[,1]*vartimestep, fsqt[,2:5]) # multiply the first column by vartimestep to give time in seconds
  fsqt <- rbind(c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"), fsqt)
  write(t(fsqt), file=paste(vardirname, "isf.txt"), ncolumns=6, sep="\t")
  
  #Do g(r) and write it out
  gr <- gr2d(tr,nbins=200,deltar=0.5,imgsize=c(varedgecutoff,varedgecutoff,varimgx-varedgecutoff,varimgy-varedgecutoff), nframes=vargofrframes)
  gr <- cbind(gr[,1], gr[,1]/varparticlesize, gr[,2]) # multiply the first column by vartimestep to give time in seconds
  gr <- rbind(c("Distance (Diameters)", "Distance (Pixels)", "g(r)"), gr)
  write(t(gr),file=paste(vardirname, "gr.txt"),ncolumns=3,sep="\t")
  
}