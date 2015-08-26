# peterpinning.r
# Peter Crowther
# Sept 2015

# Calculate number of pinned particles in each frame?
# Uses output from pintagoct.r

peterpinning = function(input){

  nparticles <- max(input[,7])
  
  cat("Progress of particle pin determination\n")
  objprogress <- txtProgressBar(min=0, max=nparticles-1, style=3)
  
  input <- input[order(input[,7]), ]   # Sorting is super fast - make use of this to avoid use of which()
  
  displacement <- matrix(data = 0, ncol=1,nrow=nparticles)
  linecounter <- 1
  for(i in 1:(nparticles-1)){
    
    if (i %% 100 == 0) {
      setTxtProgressBar(objprogress, i)
    }
    
    # Section to seperate out each frame of particle i into a matrix called thisparticle
    thisparticle <- matrix(data=0, nrow=1, ncol=8)
    breaker <- FALSE
    
    while (breaker == FALSE) {
      if (input[linecounter,7] == i){
        newrow <- input[linecounter,]
        thisparticle <- rbind(thisparticle,newrow)
        linecounter <- linecounter + 1
      }
      else{
        breaker <- TRUE
      }
    }

    
    #Get displacement of the particle in each frame
	  for (j in 2:(nrow(thisparticle)-1)){    #Skip the first row, it is blank
		  displacement[i] = displacement[i] +(((thisparticle[j, 1] - thisparticle[j+1, 1])^2 + (thisparticle[j, 2] - thisparticle[j+1, 2])^2)^(1/2))/(thisparticle[j+1, 6]-thisparticle[j, 6])
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