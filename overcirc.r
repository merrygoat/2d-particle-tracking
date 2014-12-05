# overcirc.r
# Overlays feature co-ordinates on an image
# Replaces overlayfeatures.r

overcirc = function(img,coords,rad=2){
  
  nparticles <- nrow(coords)
  
  
  
  for (i in 1:nparticles){
    img <- drawCircle(img,coords[i,1],coords[i,2],radius=rad,col=1)
  }
  
  display(img,method="raster")
  
  #return(img)
  
}