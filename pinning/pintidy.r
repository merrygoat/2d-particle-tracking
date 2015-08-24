# pintidy.r
# Ian Williams
# May 2015

# Tidies up the pintag data by removing mistaken pin identifcations
# Basically, if a particle is pinned but subsequently unpinned then it is not pinned, I have just identified it incorrectly.

pintidy = function(input){
  
  # Column 9 is the pin identification column
  # Which particles are pinned ever?
  pinids <- unique(input[input[,9]==1,7])
  # pinids holds the particle IDs for all the particles that are ever identified as pinned in the video
  
  npinned <- length(pinids)
  
  # time for some run length encoding
  # Loop on particle IDs
  for(i in 1:npinned){
    # Get this particle
    
    cat("Checking pin",i,"of",npinned,"\n")
    
    thisparticle <- input[input[,7]==pinids[i],]
    
    runlen <- rle(thisparticle[,9])
    
    # If there are more than 2 runs then the particle is not pinned. 
    # i.e. all that is allowed is the 0->1 transition. A transition back to 0 means I have messed up the identification
    
    if((length(runlen$values)>2)&(runlen$values[1]==0)){
      input[input[,7]==pinids[i],9] <- 0
    }
    
    if((length(runlen$values)>1)&(runlen$values[1]==1)){
      input[input[,7]==pinids[i],9] <- 0
    }
    
  }
  
  return(input)
  
}