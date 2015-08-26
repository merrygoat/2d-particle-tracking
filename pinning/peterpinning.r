# peterpinning.r
# Peter Crowther
# Sept 2015

# Calculate number of pinned particles in each frame?
# Uses output from pintagoct.r

peterpinning = function(input, numsplits=10){

  nparticles <- max(input[,7])
  
  cat("Progress of particle pin determination\n")
  objprogress <- txtProgressBar(min=0, max=nparticles-1, style=3)
  
  input <- input[order(input[,7]), ]   # Sorting is super fast - make use of this to avoid use of which()
  
  displacement <- matrix(data = 0, ncol=1,nrow=nparticles)
  linecounter <- 1
  for(particlenum in 1:(nparticles-1)){
    
    if (particlenum %% 100 == 0) {
      setTxtProgressBar(objprogress, particlenum)
    }
    
    # Section to seperate out each frame of particlenum into a matrix called thisparticle
    thisparticle <- matrix(data=0, nrow=1, ncol=8)
      
    repeat {
      if (input[linecounter,7] == particlenum){
        newrow <- input[linecounter,]
        thisparticle <- rbind(thisparticle,newrow)
        linecounter <- linecounter + 1
      }
      else{
        break
      }
    }
    
    #Get displacement of the particle in each frame
	  for (j in 2:(nrow(thisparticle)-1)){    #Skip the first row, it is blank
		  displacement[particlenum] = displacement[particlenum] +(((thisparticle[j, 1] - thisparticle[j+1, 1])^2 + (thisparticle[j, 2] - thisparticle[j+1, 2])^2)^(1/2))
	  }  
  }
  
  hist(displacement, breaks = seq(0,1000), xlim=c(0,50))
  
  #output <- output[2:nrow(output),]
  
  # Sort by frame rather than by particle
  #output <- output[sort.list(output[,6]),]
  close(objprogress)
  #return(output)
  return(0)
}