# gensin.r
# Ian Williams
# August 2014

# Generates points on a sin wave
# Test function for learning R

gensin = function(){
  output <- matrix(ncol=2,nrow=100)
  output[,1] <- seq(from=0,to=9.9,by=(0.1))
  output[,2] <- sin(output[,1])
  return(output)
}