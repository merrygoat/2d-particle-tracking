# sdpchpintag.r
# Ian Williams
# February 2015

# Adapted sdposchunks to tag pinned or not pinned in each time chunk
# Adds a column to trsize (or tr, I guess) that is 1 for pinned and 0 for not pinned

sdpchpintag = function(input,chunkwidth,sdthresh){
  
  nparticles <- max(input[,7])
  
  # output is: trsize + additional pin tag column
  output <- matrix(nrow=1,ncol=(ncol(input)+1))
  
  for(i in 1:nparticles){
    
    cat("Particle:",i,"of",nparticles,"\n")
    thisparticle <- input[input[,7]==i,]
    pinned <- matrix(ncol=1,nrow=nrow(thisparticle))
    
    if(nrow(thisparticle)>chunkwidth){
      #thisparticleout <- matrix(nrow=nrow(thisparticle),ncol=4)
      thisparticleposchunks <- sdposchunks(thisparticle,chunkwidth=chunkwidth)
      nchunks <- nrow(thisparticleposchunks)
      
      # Loop over the chunks
      for (j in 1:nchunks){
        # Chunk start is in col 1, chunk end is col 2, sdpos is in col 3
        if (thisparticleposchunks[j,3]<=sdthresh){
          # Particle is pinned in chunk j
          # Find all rows in this chunk interval
          w <- which((thisparticle[,6]>=thisparticleposchunks[j,1])&(thisparticle[,6]<=thisparticleposchunks[j,2]))
          pinned[w,1] <- 1
        } else{
          # Particle is not pinned in chunk j
          w <- which((thisparticle[,6]>=thisparticleposchunks[j,1])&(thisparticle[,6]<=thisparticleposchunks[j,2]))
          pinned[w,1] <- 0
        }
        
        # Need to deal with the data the exceeds the final chunk
        # Set the pinned value to the value it had in the final complete chunk
        wexceed <- which(thisparticle[,6]>thisparticleposchunks[nrow(thisparticleposchunks),2])
        wlast <- which(thisparticle[,6]==thisparticleposchunks[nrow(thisparticleposchunks),2])
        pinexceed <- pinned[wlast,1]
        pinned[wexceed,1] <- pinexceed
        
      }
      
    } else{
      # If track is less than the chunkwidth then particle cannot be pinned
      pinned[,1] <- 0
      # 0 means not pinned
    }
    
    thisparticleout <- cbind(thisparticle,pinned)
    output <- rbind(output,thisparticleout)
    
  }
  
  output <- output[2:nrow(output),]
  
  # Sort by frame rather than by particle
  output <- output[sort.list(output[,6]),]
  
  return(output)
  
}