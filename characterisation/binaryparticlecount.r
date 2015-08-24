# binaryparticlecount.r
# Ian Williams
# Sept 2014

binaryparticlecount = function(tracks){
  
  nframes <- max(tracks[,6])
  
  output <- matrix(ncol=4,nrow=nframes)
  
  # output has: frame number, total particles, big particles, small particles
  
  for (i in 1:nframes){
    
    cat("Frame:",i,"\n")
    
    output[i,1] <- i
    
    w <- which(tracks[,6]==i)
    nowtracks <- tracks[w,]
    
    output[i,2] <- nrow(nowtracks)
    output[i,3] <- sum(nowtracks[,8]==1)
    output[i,4] <- sum(nowtracks[,8]==0)
    
    
  }
  
  return(output)
  
  
}