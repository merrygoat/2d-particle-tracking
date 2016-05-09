#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

threeDconfocaltrackroutine = function(remove_drift = TRUE){
  
  # source all scripts in the current directory
  #setwd("/Users/pc9836/Documents/git/2d-particle-tracking")
  setwd("C:\\Users\\Peter\\Documents\\Uni\\PhD\\2d-particle-tracking")
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  #varfilename <- "/Volumes/WIN_DATA/Confocal/STED/15-12-21/images/FITC 19"
  varfilename <- "E:\\Confocal\\STED\\15-12-21\\processed images\\FITC 19"
  #put slash on end of dirname
  #vardirname <- "/Volumes/WIN_DATA/Confocal/STED/15-12-21/"
  vardirname <- "E:\\Confocal\\STED\\15-12-21\\"
  
  #Pretrack variables
  varimages <- 30         #How many image to read from varfilename
  varzdepth <- 2 #18         #How many z planes there are
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
  vartimestep = 60        #Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  vargofrframes = 64      #How many frames of data to analyse for the g(r)
  
  
  ### Main ###
  
  testimg <- readImage(paste(varfilename,"_t0001_z0001.png",sep=""))
  
  varimgx <- dim(testimg)[1]          #Width of the image in pixels
  varimgy <- dim(testimg)[2]          #Height of the image in pixels
  
  allparticlecount <- matrix(seq(from=1, to=varimages-1))
  allmsq  <- seq(from=1, to=varimages-1)
  allmsq <- cbind(matrix(allmsq), allmsq*vartimestep)
  allsamples <- seq(from=1, to=varimages-1)
  allsamples <- cbind(matrix(allsamples), allsamples*vartimestep)
  allfsqt <- seq(from=1, to=varimages-1)
  allfsqt <- cbind(matrix(allfsqt), allfsqt*vartimestep)
    
  for(i in vardisplacement:(varzdepth+vardisplacement)){
    cat("Frame", i, " of ", varzdepth, "\n")
    pt <- pretrack(filename=varfilename,images=varimages,diameter=vardiameter,filter=varfilter,bgavg=varbgavg, masscut=varmasscut,minimum=varminimum,chan="grey",zstack=TRUE,znumber=i)
    ptfilt <- which(pt[,1] > varedgecutoff & pt[,1] < (varimgx-varedgecutoff) & pt[,2] > varedgecutoff & pt[,2] < (varimgy-varedgecutoff))
  
    
    #Get particle count for each image
    particlecount <- rep(0, varimages-1)
    for(i in 0:varimages-1) {particlecount[i] <- sum(pt[ptfilt,6] == i)}
    allparticlecount <- cbind(allparticlecount, matrix(particlecount))
    
    tr <- iantrack(pretrack=pt[ptfilt,],maxdisp=varmaxdisp,imgsize=c(varimgx,varimgy),goodenough=10)
    if(remove_drift == TRUE) {
      # write to new array - would paste over "tr" but had trouble with "tr" not taking new values
      trnodrift <- driftremoval(tr)
    }
  
    #write(t(tr),file=paste(vardirname, "raw_coords.txt", sep=""),ncolumns=7,sep="\t")
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
  allsamples <- cbind(allsamples, rowSums(allsamples[,3:(varzdepth+2)]))
  weightedfsqt <- allfsqt
  weightedmsq <- allmsq
  
  for(i in 3:(varzdepth+2)){
    weightedfsqt[,i] <- weightedfsqt[,i]*(allsamples[,i]/allsamples[,(varzdepth+3)])
    weightedmsq[,i] <- weightedmsq[,i]*(allsamples[,i]/allsamples[,(varzdepth+3)])
  }
  weightedfsqt <- cbind(weightedfsqt, rowSums(weightedfsqt[,3:(varzdepth+2)]))
  weightedmsq <- cbind(weightedmsq, rowSums(weightedmsq[,3:(varzdepth+2)]))
  
  write(t(allparticlecount), file=paste(vardirname, "particlecount.txt", sep=""), ncolumns=1+varzdepth, sep="\t")
#  write(t(allfsqt), file=paste(vardirname, "isf.txt", sep=""), ncolumns=2+varzdepth, sep="\t")
  write(t(weightedfsqt[,c(1,2,(varzdepth+3))]), file=paste(vardirname, "weightedisf.txt", sep=""), ncolumns=3, sep="\t")
  write(t(weightedmsq[,c(1,2,(varzdepth+3))]), file=paste(vardirname, "weightedmsd.txt", sep=""),ncolumns=3,sep="\t")
  write(t(allsamples),file=paste(vardirname, "samples.txt", sep=""),ncolumns=varzdepth+3,sep="\t")

}
