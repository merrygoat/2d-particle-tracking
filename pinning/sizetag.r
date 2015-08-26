# sizetag.r
# Ian Williams
# Sept 2014
# Peter Crowther
# August 2015

# Adds column to tracked data to indicate big or small in binary sample
# 1 = big, 0 = small
# Looks at each particle in turn and uses average of brightness (col 3) to discern size
# find appropriate cut-off brightness using hist
# Modified August 2015 to remove which statment. Should now scale as O(n) as opposed to O(n^2).

sizetag = function(input,cutoff){
  
  # How many particles?
  nparticles <- max(input[,7])
  
  cat("Differentiating between large and small particles...\n")

  input <- input[order(input[,7]), ]   # Sorting is super fast - make use of this to avoid use of which()
  
  output <- matrix(nrow=nrow(input),ncol=8)
  output[,1:7] <- input
    
  brightness <- matrix(data=0, nrow=nparticles, ncol=3)     # Coulmns are summed brightness, num samples and mean

  # loop through once collecting and summing the brightnesses
  for (i in 1:nrow(input)){
    particle <- input[i,7]
    brightness[particle,1] = brightness[particle,1] + input[i,3]
    brightness[particle,2] = brightness[particle,2] + 1
  }
  
  brightness[,3] = brightness[,1]/brightness[,2]    # calculate the mean brightness for each particle
  
  # loop through again, this time assigning large and small to the particles as we go
  for (i in 1:nrow(input)){
    particle <- input[i,7]
    if (brightness[particle,3] >= cutoff){
      output[i,8] <- 1
    } else{
      output[i,8] <- 0
    }
  }
  
  return(output)
  
}