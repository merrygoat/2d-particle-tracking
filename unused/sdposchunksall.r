# sdposchunksall.r

# Applies sdposchunks to all particles in turn

sdposchunksall = function(input,chunkwidth){
  
  #nparticles <- max(input[,7])
  
  # Tabulate the particle IDs and find only those particles that hang around for at leas chunkwidth
  len <- tabulate(input[,7])
  w <- which(len > chunkwidth)
  ids <- seq(from=1, to=max(input[,7]), by=1)
  ids <- ids[w]
  
  nparticles <- length(ids)
  
  # output is: start frame, end frame, sdposition, particle number
  output <- matrix(nrow=1,ncol=4)
  
  for(i in 1:nparticles){
    
    particlenum <- ids[i]
    
    cat("Particle:",i,"of",nparticles,"\n")
    thisparticle <- input[input[,7]==particlenum,]
   
    #thisparticleout <- matrix(nrow=nrow(thisparticle),ncol=4)
    thisparticleposchunks <- sdposchunks(thisparticle,chunkwidth=chunkwidth)
    particleid <- matrix(ncol=1,nrow=nrow(thisparticleposchunks))
    particleid[,1] <- i
    thisparticleout <- cbind(thisparticleposchunks,particleid)
    output <- rbind(output,thisparticleout)
    
    
  }
  
  output <- output[2:nrow(output),]
  
}