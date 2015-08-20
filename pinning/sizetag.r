# sizetag.r
# Ian Williams
# Sept 2014

# Adds column to tracked data to indicate big or small in binary sample
# 1 = big, 0 = small
# Looks at each particle in turn and uses average of brightness (col 3) to discern size
# find appropriate cut-off brightness using hist

sizetag = function(input,cutoff){
  
  # How many particles?
  nparticles <- max(input[,7])
  
  output <- matrix(nrow=nrow(input),ncol=8)
  output[,1:7] <- input
  
  # Loop over particles
  for (i in 1:nparticles){
    cat("Particle:",i,"\n")
    w <- which(input[,7]==i)
    avgbright <- mean(input[w,3])
    if (avgbright >= cutoff){
      output[w,8] <- 1
    } else{
      output[w,8] <- 0
    }
  }
  
  return(output)
  
}