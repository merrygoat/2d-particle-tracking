# sdposchunks.r
# Ian Williams
# February 2015

# Takes a particle trajectory, splits into time-chunks of a specified duration
# Calculates the standard deviation in position for each of these chunks
# Maybe helpful in pin identification

sdposchunks = function(input,chunkwidth){
  
  startframe <- min(input[,6])
  endframe <- max(input[,6])
  nframes <- endframe - startframe + 1
  
  nchunks <- floor(nframes/chunkwidth) - 1
  
  # output holds chunk start frame, chunk end frame, sdpos in chunk
  output <- matrix(nrow=(nchunks+1),ncol=3)
  
  output[,1] <- seq(from=startframe,to=(startframe + nchunks*chunkwidth),by=chunkwidth)
  output[,2] <- seq(from=(startframe+chunkwidth-1),to=(startframe+ (nchunks+1)*chunkwidth - 1),by=chunkwidth)
  
  for (i in 1:(nchunks+1)){
    startframe <- output[i,1]
    endframe <- output[i,2]
    thischunk <- input[(input[,6]>=startframe)&(input[,6]<=endframe),]
    
    meanx <- mean(thischunk[,1])
    meany <- mean(thischunk[,2])
    
    xfrommean <- thischunk[,1] - meanx
    yfrommean <- thischunk[,2] - meany
    
    rfrommean <- sqrt(xfrommean^2 + yfrommean^2)
    
    output[i,3] <- sd(rfrommean)
  }
  
  return(output)
  
}