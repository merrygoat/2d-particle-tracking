# newvoronoi.r
# Ian Williams
# November 2014


# Does voronoi from particle co-ords for single frame
# requires library(tripack)



newvoronoi = function(input){
  
  # Load tripack
  library(tripack)
  
  # Right, let's first do the voronoi
  vor <- voronoi.mosaic(input[,1],input[,2])
  
  # And extract the polygons
  polygons <- cells(vor)
  
  # Output adds 2 columns containing number of neighbours and area
  output <- matrix(ncol=(ncol(input)+2),nrow=nrow(input))
  output[,1:ncol(input)] <- input
  
  for (i in 1:nrow(input)){
    #cat("Particle:",i,"\n")
    #output[i,(ncol(input)+1):(ncol(input)+2)] <- polygons[[i]]$center
    output[i,(ncol(input)+1)] <- length(polygons[[i]]$neighbours)
    output[i,(ncol(input)+2)] <- polygons[[i]]$area
  }
  
  return(output)
  
}