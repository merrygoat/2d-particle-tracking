# isf.r
# Computes self ISF

isf = function(tracks,interval=1,length,maxtime=0,drift=FALSE){
  
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
  # dt, real, imaginary, modulus, N
  output <- matrix(nrow=ntimes,ncol=5)
  output[,1] <- dts
  
  cat("Progress of isf\n")
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
    
    if (drift==TRUE){
      driftx <- mean(dx)
      drifty <- mean(dy)
      dx2 <- (dx - driftx)^2
      dy2 <- (dy - drifty)^2
    } else{
      dx2 <- dx^2
      dy2 <- dy^2
    }
    
    dx2 <- dx^2
    dy2 <- dy^2
    dr2 <- dx2 + dy2
    dr <- sqrt(dr2)
    cosdr <- cos((2*pi/length)*dr)
    sindr <- sin((2*pi/length)*dr)
    modulus <-  sum(besselJ((2*pi/length)*dr,0))/length(w)
  # Francesco -> Ian's version, suspected to be wrong:
  # sqrt(sum(cosdr)^2 + sum(sindr)^2)/length(w)
    
    output[i,2] <- mean(cosdr)
    output[i,3] <- mean(sindr)
    output[i,4] <- modulus
    output[i,5] <- length(w)
    
  }
  
  close(objprogress)
  
  return(output)
  
  
}