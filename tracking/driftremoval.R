#!/usr/bin/env Rscript

# Drift removal for 2D confocal images
# Peter Crowther - November 2015

driftremoval = function(trackarray, display_graphs=FALSE){

  numframes = max(trackarray[,6])
  numparticles = max(trackarray[,7])
  numobservations = length(trackarray[, 7])
  totaldisplacement = matrix(nrow=numframes, ncol=5, data=0)
  
  cat("Progress of drift removal\n")
  objprogress <- txtProgressBar(min=0, max=numparticles, style=3)

  # This section is to avoid use of a which statement which scales very poorly with high N
  # Order the data by particle number and then by frame number
  trackarray <- trackarray[order(trackarray[,7],trackarray[,6]), ]
  
  dictionary = matrix(nrow=numparticles, ncol=1, data=0)
  
  for(i in 1:(numobservations-1)) {
    if(trackarray[i, 7] != trackarray[i+1, 7]) {
      dictionary[trackarray[i+1, 7]] = i 
    }
  }
  
  for(particlenumber in 1:(numparticles-1)) {
    
    if(particlenumber%%100==0){setTxtProgressBar(objprogress, particlenumber)}
    
    particlexdisplacement = 0
    particleydisplacement = 0

    particletrack = trackarray[(dictionary[particlenumber]+1):dictionary[particlenumber+1],]
    
    if(length(particletrack) > 7) {    #make sure that we have more than one particle
      for(i in 1:(length(particletrack[,6])-1)) {
        framenumber = particletrack[i,6]
        totaldisplacement[framenumber+1,1:2] = totaldisplacement[framenumber+1,1:2] + (particletrack[i,1:2]-particletrack[i+1,1:2])
        totaldisplacement[framenumber+1,3] = totaldisplacement[framenumber+1,3] + 1
      }
    }
  }
  
  ######## Normalise displacement and Calculate cumulative displacement ############

  totaldisplacement[,1:2] = totaldisplacement[,1:2] / totaldisplacement[,3]
  totaldisplacement[1,4:5] = totaldisplacement[1,1:2]
  if (length(totaldisplacement[,1]) > 1) {
    for(i in 2:length(totaldisplacement[,1])){
      totaldisplacement[i,4:5] = totaldisplacement[i,1:2] + totaldisplacement[i-1,4:5]
    }
  }
  
  ################  and plot it #################

  if(display_graphs==TRUE){
    if(.Platform$OS.type == "unix") {x11()} else{windows()}
    plot(seq(from = 1, to = numframes), totaldisplacement[1:length(totaldisplacement[,4]),4], ylab = "cumulative x displacement")
    lines(seq(from = 1, to = numframes), totaldisplacement[1:length(totaldisplacement[,4]),4])
    if(.Platform$OS.type == "unix") {x11()} else{windows()}
    plot(seq(from = 1, to = numframes), totaldisplacement[1:length(totaldisplacement[,4]),5], ylab = "cumulative y displacement")
    lines(seq(from = 1, to = numframes), totaldisplacement[1:length(totaldisplacement[,4]),5])
  }
  
  ################ Remove drift from tracks and return array ######################
  
  for(i in 1:numframes){
    framemask = which(trackarray[,6] == i)
    for(j in 1:length(framemask)){
      trackarray[framemask[j], 1:2] = trackarray[framemask[j], 1:2] + totaldisplacement[i,4:5]
    }
  }
  
  return(trackarray)
}