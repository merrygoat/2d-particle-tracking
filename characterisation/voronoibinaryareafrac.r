# voronoibinaryareafrac.r
# Ian Williams
# Jan 2015

# Works like particle count except on voronoi tracks where edge particles have been removed
# Calculates area fraction based on total Voronoi cell area for non-edge particles
# Is hopefully a better measure of phi than a straight up count


voronoibinaryareafrac = function(tracks,sigmabig=24.3,sigmasmall=14.9){
  
  nframes <- max(tracks[,6])
  
  # Particle areas
  areabig <- pi*(sigmabig/2.0)^2
  areasmall <- pi*(sigmasmall/2.0)^2
  
  output <- matrix(ncol=8,nrow=nframes)
  
  # output has: frame number, total particles, big particles, small particles, total voronoi cell area, big area frac, small area frac total area frac
  
  for (i in 1:nframes){
    
    cat("Frame:",i,"\n")
    
    output[i,1] <- i
    
    w <- which(tracks[,6]==i)
    nowtracks <- tracks[w,]
    
    output[i,2] <- nrow(nowtracks)
    output[i,3] <- sum(nowtracks[,8]==1)
    output[i,4] <- sum(nowtracks[,8]==0)
    output[i,5] <- sum(nowtracks[,10],na.rm=TRUE)
    
    # N_big * A_big, N_small * A_small
    totbigarea <- output[i,3]*areabig
    totsmallarea <- output[i,4]*areasmall
    
    output[i,6] <- totbigarea/output[i,5]
    output[i,7] <- totsmallarea/output[i,5]
    output[i,8] <- (totbigarea + totsmallarea)/output[i,5]
    
  }
  
  return(output)
  
  
}