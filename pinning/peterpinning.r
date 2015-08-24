# peterpinning.r
# Peter Crowther
# Sept 2015

# Calculate number of pinned particles in each frame?
# Uses output from pintagoct.r

peterpinning = function(input,chunkwidth,sdthresh){

  nparticles <- max(input[,7])
  
  # output is: trsize + additional pin tag column
  output <- matrix(nrow=1,ncol=(ncol(input)+1))
  
  cat("Progress of particle pin determination\n")
  objprogress <- txtProgressBar(min=0, max=nparticles-1, style=3)
  
  displacement <- matrix(ncol=1,nrow=nrow(nparticles))
  
  for(i in 1:nparticles){
    
    if (i %% 10 == 0) {
      setTxtProgressBar(objprogress, i)
    }
	
    thisparticle <- input[input[,7]==i,]
    pinned <- matrix(ncol=1,nrow=nrow(thisparticle))
    	
	#Get displacement of the particle in each frame
	for (j in 1:nrow(thisparticle)){
		displacement[i] = displacement[i] +(((thisparticle[j][1] - thisparticle[j+1][1])^2 + (thisparticle[j][1] - thisparticle[j+1][1])^2)^(1/2))/(thisparticle[j][6]-thisparticle[j][6])
	}  
  }
  
  hist(displacement)
  
  #output <- output[2:nrow(output),]
  
  # Sort by frame rather than by particle
  #output <- output[sort.list(output[,6]),]
  close(objprogress)
  #return(output)
  return(0)
}