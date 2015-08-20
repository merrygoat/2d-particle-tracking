# pinsintime.r
# Ian Williams
# Sept 2014

# Calculate number of pinned particles in each frame?
# Uses output from pintagoct.r

pinsintime = function(tracks){
  
  nframes <- max(tracks[,6])
  
  output <- matrix(ncol=7,nrow=nframes)

  # output has: frame number, total particles, big particles, small particles, pinned particles, big pinned, small pinned
  
  for (i in 1:nframes){
    
    cat("Frame:",i,"\n")
    
    output[i,1] <- i
    
    w <- which(tracks[,6]==i)
    nowtracks <- tracks[w,]
    
    output[i,2] <- nrow(nowtracks)
    output[i,3] <- sum(nowtracks[,8]==1)
    output[i,4] <- sum(nowtracks[,8]==0)
    output[i,5] <- sum(nowtracks[,9]==1)
    output[i,6] <- sum(nowtracks[,9]==1&nowtracks[,8]==1)
    output[i,7] <- sum(nowtracks[,9]==1&nowtracks[,8]==0)
    
    
  }
  
  return(output)
  
  
}