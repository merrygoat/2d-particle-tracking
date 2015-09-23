# pretrack.r
# Reads in, processes and pretracks a stack of tiffs numbered sequentially
# Require lowpass.r and feature.r
# Ian Williams
# 23/09/2013
#
# Requires library(EBImage)
#
# There are many parameters here:
#    filename: Path and filename stub string, with double backslashes for R
#    images: number of images in time series
#    crop: vector of 4 numbers indicating cropping of form c(xmin,xmax,ymin,ymax)
#          defaults to c(0,0,0,0) which indicates no cropping
#    rotang: angle for rotation in degrees. Rotation applied in anti-clockwise direction and occurs before cropping
#          defaults to zero
#    bgavg: size of background averaging mask in image processing, must be odd integer
#          defaults to 5
#    filter: range of the lowpass filter (called lobject in lowpass.r)
#    diameter: approximate diameter of particles (slight overestimate is best)
#    separation: minimum separation for particles
#          if two features are within separation of each other the brightest is retained
#    masscut: minimum total brightness for good features (much less than JC code)
#    minimum: minimum peak brightness for good features (value between 0 and 1)
#          feature only considered if brightest pixel is at least this bright
#
#
# Output is 5 column matrix with one row per particle per frame
# Can be written to text with
# write(t(f),file="C:\\Users\\chxiw\\Dropbox\\FoolingR\\r_pt.dat",ncolumns=4,sep="\t")
#
# Output is:
# x position, y position, brightness, radius of gyration, frame number
#
# Note: Eccentricity is not computed because of laziness and the fact that I rarely use it for anything

pretrack = function(filename,images,crop=c(0,0,0,0),rotang=0,filter,bgavg=5,diameter,separation=0,masscut,minimum,chan="grey"){
  
  # In contrast with old IDL code, I will read in and handle each image separately in order
  # rather than reading in whole sequence prior to processing
  
  # Order of operations is:
  # 1. Read image
  # 2. Crop image
  # 3. Lowpass image
  # 4. Find features
  # 5. Append info to pretrack array
  
  # Assumes images is the total number of images, first image numbered 000000
  # So we loop up to images-1
  
  # Make output array of 5 columns to store:
  # x, y, mass, r^2, frame
  # 1 row for now, but we will ignore the first row at the end
  output <- matrix(nrow=1,ncol=6)
  #cat(output) #currently says NA NA NA NA NA
    
  cat("Progress of particle identification\n")

  objprogress <- txtProgressBar(min=0, max=images-1, style=3)
    
  for (i in 0:(images-1)){
    setTxtProgressBar(objprogress, i)  
    

    # Build the filepath and name for this image
    #cat("Reading image ",i,"\n")
    thisimagename <- paste(filename,formatC(i,flag="0",digits=3),".png",sep="")
    #cat(thisimagename,"\n")
    thisimage <- channel(readImage(thisimagename), chan)
    
    
    # Rotate this image if you so desire
    if (rotang!=0){
      thisimage <- rotate(thisimage,rotang)
    }
        
    # Crop this image if you so desire
    # By default crop parameters are set to zero, meaning no crop
    if (sum(crop) != 0){
      thisimage <- thisimage[crop[1]:crop[2],crop[3]:crop[4]]      
    }
    
    # Lowpass this frame
    lp <- lowpass(image=thisimage,lobject=filter,bgavg=bgavg)
    
    # Find features
    f <- feature(image=lp,diameter=diameter,separation=separation,minimum=minimum,masscut=masscut)
    
    # Add another column to f to hold the frame number
    # How many particles in this frame?
    thislength <- nrow(f)
    framenum <- matrix(ncol=1,nrow=thislength)
    framenum <- i
    f <- cbind(f,framenum,deparse.level=0)
    
    # Append this to the output matrix
    output <- rbind(output,f)
    
  }
  
  close(objprogress)
  
  # How long is output?
  length <- nrow(output)
  # Remove first row because it is NA NA NA NA NA
  output <- output[2:length,]
  
  return(output)
  
  
}