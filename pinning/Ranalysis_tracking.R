#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D brightfield images of pinned particles

pinningtrackroutine = function(ptdone=FALSE) {
  
  # Setup r scripts
  setwd("/Users/pc9836/Documents/git/2d-particle-tracking")
  filelist <- c("pre_tracking/lowpass.r", "pre_tracking/feature.r", "pre_tracking/pretrack.r", "tracking/iantrack.r", 
                "pinning/sizetag.r", "pinning/sdpchpintag.r", "pinning/sdposchunks.r", #"pinning/pintidy.r",
                "pinning/pinsintime.r", #"characterisation/binaryparticlecount.r",
                "characterisation/gr2d.r", "characterisation/msd.r", "characterisation/isf.r",
                "characterisation/newvoroseries.r", "characterisation/voronoibinaryareafrac.r", "characterisation/neighboursbinary.r" )
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)
  
  #File directory variables
  istiffstack <- FALSE     #Set to true if the images are read in as a single compound tiff image
  varfilename <- "/Volumes/WIN_DATA/Pinning/Slide_12/10_mins/sample12-2pm_0000"
  vardirname <- "/Volumes/WIN_DATA/Pinning/Slide_12/10_mins/"
  setwd(vardirname)
  
  #Setting variables
  varcutoff <- 26      # the diameter cutoff between large and small particles
  
  #Pretrack variables
  varimages <- 100        #How many image to read from varfilename
  vardiameter <-15        #Particle diamter - used in particle identification
  varfilter <-11          #Parameter for lowpass filter
  varbgavg <- 13          #Parameter for lowpass filter
  varmasscut <- 8         #Lowest integrated brightness for particle
  varminimum <- 0.2       #Lowest pixel value for center of particle
  
  #Track variables
  varedgecutoff <- 15     #Cuts off this many pixels from each edge of the image in all data output - this is because particle identification is bad around the edges.
  varmaxdisp <- 5         #Used in tracking - the maximum allowed interframe displacement
  
  #Pinning varibales
  varchunkwidth <- 1000   # The time window used to look for pinning. Alter as a function of the relaxation time of the system.
  versdthresh <- 0.3      # The threshold below which things are pinned
  
  # Science variables
  varbigparticlesize = 24.3    #Used as the wavevector for isf
  varsmallparticlesize = 14.9    #Used as the wavevector for isf
  
  ########################### Main ###########################
  
  if (istiffstack == TRUE) {
    cat("Reading image file.\n")
    imgtiffstack <- suppressWarnings(readImage(paste(varfilename,".tif",sep="")))   #Supress to avoid warnings about unreadable metadata
  }
  else {
    imgtiffstack <- suppressWarnings(readImage(paste(varfilename,"0000.tif",sep=""))) #If no tiff stack, read in one image anyway to get the dimensions.
  }
  
  varimgx <- dim(imgtiffstack)[1]          #Width of the image in pixels
  varimgy <- dim(imgtiffstack)[2]          #Height of the image in pixels
  
  # Particle identification
  if (ptdone == FALSE) {
    pt <- pretrack(filename=varfilename,images=varimages,diameter=vardiameter,masscut=varmasscut,minimum=varminimum,filter=varfilter,bgavg=varbgavg)
    pt[,6] <- pt[,6] - 1
    write(t(pt),file="pretrack.dat",ncolumns=6,sep="\t")
  }
  else {
    pt <- read.table("pretrack.dat",sep="\t")
    pt <- data.matrix(pt)
  }
  
  # Filter the data to cut out particles which are withing varedgecutoff pixels of the edge of the image.
  ptfilt <- which(pt[,1] > varedgecutoff & pt[,1] < (varimgx-varedgecutoff) & pt[,2] > varedgecutoff & pt[,2] < (varimgy-varedgecutoff))

  # Particle tracking
  tr <- iantrack(pretrack=pt[ptfilt,],maxdisp=varmaxdisp,imgsize=c(varimgx,varimgy),goodenough=2)
  
  # Cutoff between big and small particles. Do a histo here to determine the cutoff.
  hist(tr[,3],breaks = seq(0,1000), xlim=c(10,35))
  trsize <- sizetag(tr,cutoff=varcutoff)
  write(t(trsize),file="track.dat",ncolumns=8,sep="\t")
  # 1 = big, 0 = small
  
  # Determine what is pinned
  pntag <- sdpchpintag(trsize,chunkwidth=varchunkwidth,sdthresh=varsdthresh)
  
  pntag <- pintidy(pntag)
  write(t(pntag),file="track_sizepintag.dat",ncolumns=9,sep="\t")
  
  pntime <- pinsintime(pntag)
  write(t(pntime),file="pinsintime.dat",ncolumns=7,sep="\t")
  
  particles <- binaryparticlecount(trsize)
  write(t(particles),file="particlecount.dat",ncolumns=4,sep="\t")
  
  # Split tracks up a bit
  #######################
  
  w <- which(pntag[,8]==1)
  trbig <- pntag[w,]
  #write(t(trbig),file="track_big.dat",ncolumns=8,sep="\t")
  
  w <- which(pntag[,8]==0)
  trsmall <- pntag[w,]
  #write(t(trsmall),file="track_small.dat",ncolumns=8,sep="\t")
  
  w <- which(pntag[,9]==1)
  trpins <- pntag[w,]
  #write(t(trpins),file="track_pins.dat",ncolumns=9,sep="\t")
  
  w <- which(pntag[,9]==0)
  trnonpin <- pntag[w,]
  #write(t(trnonpin),file="track_unpinned.dat",ncolumns=9,sep="\t")
  
  w <- which(trbig[,9]==0)
  trbignonpin <- trbig[w,]
  #write(t(trbignonpin),file="track_unpinned_big.dat",ncolumns=9,sep="\t")
  
  w <- which(trsmall[,9]==0)
  trsmallnonpin <- trsmall[w,]
  #write(t(trsmallnonpin),file="track_unpinned_small.dat",ncolumns=9,sep="\t")
  
  # Some static and dynamic analysis
  ###################################
  w <- which(trsize[,6]<5)
  gr <- gr2d(trsize[w,],nbins=500,deltar=0.1,imgsize=c(varedgecutoff,varedgecutoff,varimgx-varedgecutoff,varimgy-varedgecutoff),binary=TRUE)
  write(t(gr),file="grbinary.dat",ncolumns=4,sep="\t")
  
  msdbig <- msd(trbignonpin,drift=TRUE)
  write(t(msdbig),file="msd_unpinned_big.dat",ncolumns=7,sep="\t")
  
  msdsmall <- msd(trsmallnonpin,drift=TRUE)
  write(t(msdsmall),file="msd_unpinned_small.dat",ncolumns=7,sep="\t")
  
  isfbig <- isf(trbig,length=varbigparticlesize,drift=TRUE)
  write(t(isfbig),file="isf_all_big.dat",ncolumns=5,sep="\t")
  
  isfsmall <- isf(trsmall,length=varsmallparticlesize,drift=TRUE)
  write(t(isfsmall),file="isf_all_small.dat",ncolumns=5,sep="\t")
  
  voro <- newvoroseries(trsize,imgsize=c(varedgecutoff,varedgecutoff,varimgx-varedgecutoff,varimgy-varedgecutoff),edgerm=50)
  write(t(voro),file="tr_voronoi.dat",ncolumns=10,sep="\t")
  
  vorphi <- voronoibinaryareafrac(voro)
  write(t(vorphi),file="areafrac_voronoi.dat",ncolumns=8,sep="\t")
  
  neigh <- neighboursbinary(voro)
  write(t(neigh),file="binaryneighbourscount.dat",ncolumns=15,sep="\t")
}