# psi6loc.r
# Ian Williams
# October 2013
#
# Computes local psi_6 for single frame
# Relies on tripack

psi6loc = function(xcoords,ycoords){
  
  # Load tripack
  library(tripack)
  
  # First do the triangulation
  triangles <- tri.mesh(xcoords,ycoords)
  
  # Extract neighbours
  nbours <- neighbours(triangles)
  # This is a list containing the neighbours for each particle
  
  # Output matrix holds:
  # xcoord, ycoord, psi_6 real, psi_6 imaginary, psi_6 modulus
  output <- matrix(nrow=1,ncol=5)
  
  # How many particles?
  nparticle <- length(nbours)
  
  # Loop over particles
  for (i in 1:nparticle){
    # Get the neighbours of particle i
    thisnbours <- nbours[i]
    # Unlist this
    thisnbours <- unlist(thisnbours)
    # Now this is a vector with the indices of neighbouring particles
    
    # How many neighbours do we have?
    nnbours <- length(thisnbours)
    
    # Make array to hold the distances between particle i and all its neighbours and the angles
    neighdisps <- matrix(nrow=nnbours,ncol=4)
    # Columns here are: xdisp, ydisp, r, theta
    
    # Reset cosine and sine totals
    costotal <- 0
    sintotal <- 0
    # These will hold real and imaginary parts of psi_6 calculation
    
    # Loop over neighbours
    for (j in 1:nnbours){
      # Get x and y displacement of neighbour j from particle i
      neighdisps[j,1] <- xcoords[thisnbours[j]] - xcoords[i]
      neighdisps[j,2] <- ycoords[thisnbours[j]] - ycoords[i]
      
      # Find scalar distance between particles
      neighdisps[j,3] <- sqrt(neighdisps[j,1]^2 + neighdisps[j,2]^2)
      
      # Finding angle is a little involved
      basetheta <- atan(abs(neighdisps[j,2]/neighdisps[j,1]))
      theta <- 0.0
      
      # Checking quadrants
      if ((neighdisps[j,1]==0)&(neighdisps[j,2]>0)){theta <- pi/2}
      if ((neighdisps[j,1]==0)&(neighdisps[j,2]<0)){theta <- (3*pi)/2}
      if ((neighdisps[j,1]>0)&(neighdisps[j,2]==0)){theta <- 0}
      if ((neighdisps[j,1]<0)&(neighdisps[j,2]==0)){theta <- pi}
      if ((neighdisps[j,1]>0)&(neighdisps[j,2]>0)){theta <- basetheta}
      if ((neighdisps[j,1]<0)&(neighdisps[j,2]>0)){theta <- pi - basetheta}
      if ((neighdisps[j,1]<0)&(neighdisps[j,2]<0)){theta <- basetheta + pi}
      if ((neighdisps[j,1]>0)&(neighdisps[j,2]<0)){theta <- (2*pi) - basetheta}
      
      neighdisps[j,4] <- theta
      
      # Add to sin and cos totals
      costotal <- costotal + cos(6*theta)
      sintotal <- sintotal + sin(6*theta)
      
    }
    
    # Contribution to output
    thisoutdata <- matrix(nrow=1,ncol=5)
    
    # Coords in first two columns
    thisoutdata[,1] <- xcoords[i]
    thisoutdata[,2] <- ycoords[i]
    
    thisoutdata[,3] <- costotal/nnbours
    thisoutdata[,4] <- sintotal/nnbours
    
    thisoutdata[,5] <- sqrt(thisoutdata[,3]^2 + thisoutdata[,4]^2)
    
    output <- rbind(output,thisoutdata)
    
  }
  
  # Remove first row of output as it is NA
  length <- nrow(output)
  # Remove first row because it is NA NA NA NA NA
  output <- output[2:length,]
  
  return(output)
  
}