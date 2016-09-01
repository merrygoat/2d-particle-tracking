#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

threeDconfocaltrackroutine = function(remove_drift = TRUE){
  
  # source all scripts in the current directory
  setwd("/Users/pc9836/Documents/git/2d-particle-tracking")
  #setwd("C:\\Users\\Peter\\Documents\\Uni\\PhD\\2d-particle-tracking")
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  varfilename <- "/Volumes/WIN_DATA/Confocal/STED/Soft\ spheres/16-08-25/v/images/v_"
  #varfilename <- "E:\\Confocal\\STED\\15-12-21\\processed images\\FITC 19"
  #put slash on end of dirname
  vardirname <- "/Volumes/WIN_DATA/Confocal/STED/Soft\ spheres/16-08-25/v/"
  #vardirname <- "E:\\Confocal\\STED\\15-12-21\\"
  
  #Pretrack variables
  varimages <- 250         #How many image to read from varfilename
  varzdepth <- 22         #How many z planes there are
  vardiameter <- 11       #Particle diameter - used in particle identification
  varfilter <- 11         #Parameter for lowpass filter
  varbgavg <- 11          #Parameter for lowpass filter
  varmasscut <- 0.1         #Lowest integrated brightness for particle
  varminimum <- 0.01       #Lowest pixel value for center of particle
  
  #Track variables
  varedgecutoff <- 10     #Cuts off this many pixels from each edge of the image in all data output - this is because particle identification is bad around the edges.
  varmaxdisp <- 5         #Used in tracking - the maximum allowed interframe displacement
  goodenough <- 5        #Minimum acceptable trajectory length
  
  #Other variables that I can't think of a title for
  varparticlesize = 15.5    #Used as the wavevector for isf
  vartimestep = 1.05        #Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  vargofrframes = 2      #How many frames of data to analyse for the g(r)
  
  ### Variable check ###
  
  if(varimages < 2) {
    cat("Error, varimages too small. At least two timesteps must be analysed.\n")
    quit()
  }

  if(goodenough > varimages) {
    cat("Good enough is too small. Setting to varimages.\n")
    goodenough <- varimages
  }
  
  ### Main ###
  
  testimg <- readImage(paste(varfilename,"_t0001_z0001.png",sep=""))
  
  varimgx <- dim(testimg)[1]          #Width of the image in pixels
  varimgy <- dim(testimg)[2]          #Height of the image in pixels
  
  allparticlecount <- matrix(seq(from=0, to=varimages-1))
  allmsq  <- seq(from=1, to=varimages-1)
  allmsq <- cbind(matrix(allmsq), allmsq*vartimestep)
  allsamples <- seq(from=1, to=varimages-1)
  allsamples <- cbind(matrix(allsamples), allsamples*vartimestep)
  allfsqt <- seq(from=1, to=varimages-1)
  allfsqt <- cbind(matrix(allfsqt), allfsqt*vartimestep)
  
  # Allows slicing of the z-stack
  vardisplacement <- 1
  
  for(i in vardisplacement:(varzdepth+vardisplacement-1)){
    cat("Frame", i, " of ", varzdepth, "\n")
    pt <- pretrack(filename=varfilename,images=varimages,diameter=vardiameter,filter=varfilter,bgavg=varbgavg, masscut=varmasscut,minimum=varminimum,chan="grey",zstack=TRUE,znumber=i)
    ptfilt <- which(pt[,1] > varedgecutoff & pt[,1] < (varimgx-varedgecutoff) & pt[,2] > varedgecutoff & pt[,2] < (varimgy-varedgecutoff))
    
    #Get particle count for each image
    particlecount <- rep(0, varimages)
    for(j in 1:varimages) {particlecount[j] <- sum(pt[ptfilt,6] == (j-1))}
    allparticlecount <- cbind(allparticlecount, matrix(particlecount))
    
    tr <- iantrack(pretrack=pt[ptfilt,],maxdisp=varmaxdisp,imgsize=c(varimgx,varimgy),goodenough=1)
    if(remove_drift == TRUE) {
      # write to new array - would paste over "tr" but had trouble with "tr" not taking new values
      trnodrift <- driftremoval(tr)
    }
  
    write(t(tr),file=paste(vardirname, "slice_", i, "_raw_coords.txt", sep=""),ncolumns=7,sep="\t")
    #Don't need to filter tr as it was already filtered from pt
    #Do msd and write it out
    if(remove_drift == TRUE) 
      {msq <- msd(trnodrift)}
    else
      {msq <- msd(tr)}
    allmsq <- cbind(allmsq, msq[,6]/varparticlesize)  # multiply the first column by vartimestep to give time in seconds, divide the 6th column to give msd in dimeters
    allsamples <- cbind(allsamples, msq[,7])
    
    #Do isf and write it out
    if(remove_drift == TRUE) 
      {fsqt <- isf(trnodrift,length=varparticlesize)}
    else
      {fsqt <- isf(tr,length=varparticlesize)}
    allfsqt <- cbind(allfsqt, fsqt[,4]) # multiply the first column by vartimestep to give time in seconds
    #fsqt <- rbind(c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"), fsqt)
  }
  
  # Post processing of data files
  allsamples <- cbind(allsamples, rowSums(allsamples[,3:(varzdepth+2), drop=FALSE]))
  weightedfsqt <- allfsqt
  weightedmsq <- allmsq
  
  for(i in 3:(varzdepth+2)){
    weightedfsqt[,i] <- weightedfsqt[,i]*(allsamples[,i]/allsamples[,(varzdepth+3)])
    weightedmsq[,i] <- weightedmsq[,i]*(allsamples[,i]/allsamples[,(varzdepth+3)])
  }
  weightedfsqt <- cbind(weightedfsqt, rowSums(weightedfsqt[,3:(varzdepth+2), drop=FALSE]))
  weightedmsq <- cbind(weightedmsq, rowSums(weightedmsq[,3:(varzdepth+2), drop=FALSE]))
  
  write(t(allparticlecount), file=paste(vardirname, "particlecount.txt", sep=""), ncolumns=1+varzdepth, sep="\t")
  write(t(weightedfsqt[,c(1,2,(varzdepth+3))]), file=paste(vardirname, "weightedisf.txt", sep=""), ncolumns=3, sep="\t")
  write(t(weightedmsq[,c(1,2,(varzdepth+3))]), file=paste(vardirname, "weightedmsd.txt", sep=""),ncolumns=3,sep="\t")
  write(t(allsamples),file=paste(vardirname, "samples.txt", sep=""),ncolumns=varzdepth+3,sep="\t")

}
