# psi6loc.r
# Ian Williams
# October 2013
#
# Computes local psi_6 for single frame
# Relies on tripack

psi6loc = function(xcoords, ycoords){
  
  # Load tripack
  library(tripack)

  # First do the triangulation
  triangles <- tri.mesh(xcoords,ycoords)
  
  # Extract a list containing the neighbours for each particle
  neighbours <- neighbours(triangles)
  num_particles <- length(neighbours)
  
  # Output matrix holds:
  # xcoord, ycoord, psi_6 real, psi_6 imaginary, psi_6 modulus
  output <- matrix(nrow=num_particles, ncol=5)

  # Loop over particles
  for (i in 1:num_particles){
    # Get the neighbours of particle i
    i_neighbours <- unlist(neighbours[i])
    num_neighbours <- length(i_neighbours)
    
    # Make array to hold the distances between particle i and all its neighbours and the angles
    # Columns are: x_displacement, y_dispisplacemet, r, theta
    neighbour_displacements <- matrix(nrow=num_neighbours,ncol=4)
    
    # Reset cosine and sine totals
    costotal <- 0
    sintotal <- 0
    # These will hold real and imaginary parts of psi_6 calculation
    
    # Loop over neighbours
    for (j in 1:num_neighbours){
      # Get x and y displacement of neighbour j from particle i
      neighbour_displacements[j,1] <- xcoords[i_neighbours[j]] - xcoords[i]
      neighbour_displacements[j,2] <- ycoords[i_neighbours[j]] - ycoords[i]
      
      # Find scalar distance between particles
      neighbour_displacements[j,3] <- sqrt(neighbour_displacements[j,1]^2 + neighbour_displacements[j,2]^2)
      
      # Finding angle is a little involved
      basetheta <- atan(abs(neighbour_displacements[j,2]/neighbour_displacements[j,1]))
      theta <- 0.0
      
      # Checking quadrants
      if ((neighbour_displacements[j,1]==0)      & (neighbour_displacements[j,2]>0))  {theta <- pi/2}
      else if ((neighbour_displacements[j,1]==0) & (neighbour_displacements[j,2]<0))  {theta <- (3*pi)/2}
      else if ((neighbour_displacements[j,1]>0)  & (neighbour_displacements[j,2]==0)) {theta <- 0}
      else if ((neighbour_displacements[j,1]<0)  & (neighbour_displacements[j,2]==0)) {theta <- pi}
      else if ((neighbour_displacements[j,1]>0)  & (neighbour_displacements[j,2]>0))  {theta <- basetheta}
      else if ((neighbour_displacements[j,1]<0)  & (neighbour_displacements[j,2]>0))  {theta <- pi - basetheta}
      else if ((neighbour_displacements[j,1]<0)  & (neighbour_displacements[j,2]<0))  {theta <- basetheta + pi}
      else if ((neighbour_displacements[j,1]>0)  & (neighbour_displacements[j,2]<0))  {theta <- (2*pi) - basetheta}
      
      neighbour_displacements[j,4] <- theta
      
      # Add to sin and cos totals
      costotal <- costotal + cos(6*theta)
      sintotal <- sintotal + sin(6*theta)
      
    }
    
    # Coords in first two columns
    output[i, 1] <- xcoords[i]
    output[i, 2] <- ycoords[i]
    output[i, 3] <- costotal/num_neighbours
    output[i, 4] <- sintotal/num_neighbours
    output[i, 5] <- sqrt(output[i, 3]^2 + output[i, 4]^2)
    
  }

  return(output)
  
}