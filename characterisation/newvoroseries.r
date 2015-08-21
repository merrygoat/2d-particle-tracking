# newvoroseries.r
# Ian Williams
# November 2014

newvoroseries = function(input,imgsize=c(768,768),edgerm=20){
  
  nframes <- max(input[,6])
  output <- matrix(ncol=(ncol(input)+2),nrow=1)
  
  for(i in 0:nframes){
    cat("Frame:",i,"\n")
    output <- rbind(output,newvoronoi(input[input[,6]==i,]))
    
  }
  
  output<-output[2:nrow(output),]
  
  # Remove all particles within edgerm of the image edge as their cells are dumb
  w <- which((output[,1]>edgerm)&(output[,1]<(imgsize[1]-edgerm))&(output[,2]>edgerm)&(output[,2]<(imgsize[2]-edgerm)))
  output <- output[w,]
  
  return(output)
  
  
}