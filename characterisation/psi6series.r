# psi6series.r
# Ian Williams
# October 2013
#
# Applies psi6loc to all frames in a video

psi6series = function(data){
  
  # Assume first frame is 0 and we increment by 1
  # Number of frames
  nframes <- max(data[,6])
  
  # Output is 6 column matrix
  # Each row is a particle
  # Columns are: xcoord, ycoord, psi_6 real, psi_6 imaginary, psi_6 modulus, frame number
  output <- matrix(ncol=6,nrow=1)
  
  # Loop over frames
  for (i in 0:nframes){
    
    cat("Computing... Frame: ",i,"\n")
    
    # Get this frame data
    w <- which(data[,6]==i)
    thisframe <- data[w,]
    
    # Apply vorpolys
    thisframepsi6 <- psi6loc(thisframe[,1],thisframe[,2])
    
    # How many particles in this frame?
    nparticles <- nrow(thisframepsi6)
    
    # Make frame number matrix
    framenum <- matrix(nrow=nparticles,ncol=1)
    framenum[,1] <- i
    
    # Bind the two
    dataforoutput <- cbind(thisframepsi6,framenum)
    
    # Bind to output
    output <- rbind(output,dataforoutput)
    
  }
  
  # Remove first row of output as it is NA
  length <- nrow(output)
  # Remove first row because it is NA NA NA NA NA
  output <- output[2:length,]
  
  return(output)
  
}