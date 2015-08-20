# qoft.r
# Ian Williams
# May 2014

# Compute Q(t) for a given dt used in finding chi_4
# cutoff typically 0.3 sigma

qoft = function(tracks,dt,cutoff){
  
  # Max frame number
  nframes <- max(tracks[,6])
  
  # Sort the tracks to be arranged by particle ID rather than frame
  tracks <- tracks[sort.list(tracks[,7]),]
  
  # Assume number of particles is fixed
  nparticles <- length(unique(tracks[,7]))
  
  # Make output
  output <- matrix(nrow=(nframes-dt),ncol=2)
  
  # Shift tracks
  shiftedtracks <- shift(tracks,dt)
  
  # Loop over nframes-dt
  for (i in 0:(nframes-dt)){
    
    cat("Interval :",dt,"\t","Frame :",i,"\n")
    
    # Want co-ords separated by interval i
    w <- which((abs(shiftedtracks[,6]-tracks[,6])==dt)&(tracks[,6]==i))
    
    dxs2 <- (shiftedtracks[w,1] - tracks[w,1])^2
    dys2 <- (shiftedtracks[w,2] - tracks[w,2])^2
    drs <- sqrt(dxs2 + dys2)
    
    # This is all pairs of frames
    w <- which(drs <= cutoff)
    
    output[i,2] <- length(w)/nparticles
    output[i,1] <- i+dt
    
  }
  
  return(output)
  
  
}