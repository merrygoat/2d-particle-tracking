# msd.r
# computes MSD from tracked data

msd = function(tracks,interval=1,maxtime=0,drift=FALSE){
  
  # Find maximum time
  if (maxtime==0){
    tmax <- max(tracks[,6])
  } else { tmax <- maxtime }
  
  # Make output time interval list
  ntimes <- floor(tmax/interval)
  maxdt <- ntimes*interval
  dts <- seq(from=interval,to=maxdt,by=interval)
  
  # Sort the tracks to be arranged by particle ID rather than frame
  tracks <- tracks[sort.list(tracks[,7]),]
  
  # Output array holds:
  # dt, <x>, <y>, <x^2>, <y^2>, MSD, N
  output <- matrix(nrow=ntimes,ncol=7)
  output[,1] <- dts
  
  cat("Progress of msd\n")
  objprogress <- txtProgressBar(min=0, max=ntimes, style=3)
  
  # Loop over time intervals
  for (i in 1:ntimes){
    
    setTxtProgressBar(objprogress, i)  
    
    # We are considering time interval dt on this loop iteration
    dt <- output[i,1]
    
    #cat("Time interval :",dt,"\n")
    
    # Create shifted array
    shiftedtracks <- shift(tracks,dt)
    
    # Compare!
    # We want data for the same particle ID but separated by dt
    w <- which((abs(shiftedtracks[,6]-tracks[,6])==dt)&(shiftedtracks[,7]==tracks[,7]))
    
    # x and y displacements
    dx <- shiftedtracks[w,1] - tracks[w,1]
    dy <- shiftedtracks[w,2] - tracks[w,2]
    
    output[i,2] <- mean(dx)
    output[i,3] <- mean(dy)
    
    if (drift==FALSE){
      dx2 <- dx^2
      dy2 <- dy^2
      dr2 <- dx2 + dy2
    } else{
      dx2 <- (dx - output[i,2])^2
      dy2 <- (dy - output[i,3])^2
      dr2 <- dx2 + dy2
    }
    
    output[i,4] <- mean(dx2)
    output[i,5] <- mean(dy2)
    output[i,6] <- mean(dr2)
    output[i,7] <- length(w)
        
    
  }
  
  close(objprogress)
  
  return(output)
  
}