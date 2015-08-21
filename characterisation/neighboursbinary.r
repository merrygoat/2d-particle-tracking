# neighboursbinary.r
# Ian Williams
# Feb 2015

# More simple particle counting for Voronoi neighbours in binary samples
# Developed for pinning experiments
# Take Voronoi track inputs with particle size tag in column 8 and neighbours in column 9
# This only considers those particles with bounded Voronoi cells, i.e. have previously thrown out edges

neighboursbinary = function(data){
  
  # How many frames
  nframes <- max(data[,6]) + 1
  
  # Output is frame num, N_big3, N_big4, N_big5, N_big6, N_big7, N_big8, N_big9, N_small3, N_small4, N_small5, N_small6, N_small7, N_small8, N_small9
  output <- matrix(nrow=nframes,ncol=15)
  
  # Loop over frames
  for (i in 1:nframes){
    output[i,1] <- (i-1)
    
    # Get this frame data
    w <- which(data[,6]==(i-1))
    thisframedata <- data[w,]
    
    # Split into big and small
    w <- which(thisframedata[,8]==1)
    thisframebig <- thisframedata[w,]
    w <- which(thisframedata[,8]==0)
    thisframesmall <- thisframedata[w,]
    
    histbig <- hist(thisframebig[,9],breaks=seq(from=0,to=15,by=1),plot=FALSE)
    output[i,2] <- histbig$counts[3]
    output[i,3] <- histbig$counts[4]
    output[i,4] <- histbig$counts[5]
    output[i,5] <- histbig$counts[6]
    output[i,6] <- histbig$counts[7]
    output[i,7] <- histbig$counts[8]
    output[i,8] <- histbig$counts[9]
    
    histsmall <- hist(thisframesmall[,9],breaks=seq(from=0,to=15,by=1),plot=FALSE)
    output[i,9] <- histsmall$counts[3]
    output[i,10] <- histsmall$counts[4]
    output[i,11] <- histsmall$counts[5]
    output[i,12] <- histsmall$counts[6]
    output[i,13] <- histsmall$counts[7]
    output[i,14] <- histsmall$counts[8]
    output[i,15] <- histsmall$counts[9]
    
    cat("Frame:",i,"\n")
    
  }
  
  return(output)
  
}