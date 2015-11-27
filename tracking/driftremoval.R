#!/usr/bin/env Rscript

# Drift removal for 2D confocal images
# Peter Crowther - November 2015

driftremoval = function(trackarray){

  numframes = max(trackarray[,6])
  numparticles = max(trackarray[,7])
  totaldisplacement = matrix(nrow=numframes, ncol=5, data=0)
  
  for(particlenumber in 1:numparticles) {
    
    particlexdisplacement = 0
    particleydisplacement = 0
    
    particlemask = which(trackarray[,7] == particlenumber)
    particletrack = trackarray[particlemask,]
    
    for(i in 1:length(particletrack[,6])-1) {
      framenumber = particletrack[i,6]
      totaldisplacement[framenumber+1,1:2] = totaldisplacement[framenumber+1,1:2] + (particletrack[i,1:2]-particletrack[i+1,1:2])
      totaldisplacement[framenumber+1,3] = totaldisplacement[framenumber+1,3] + 1
    }
  
  }
  totaldisplacement[,1:2] = totaldisplacement[,1:2] / totaldisplacement[,3]
  totaldisplacement[1,4:5] = totaldisplacement[1,1:2]
  for(i in 2:(length(totaldisplacement[,1])-1)){
    totaldisplacement[i,4:5] = totaldisplacement[i,1:2] + totaldisplacement[i-1,4:5]
  }
  windows()
  plot(seq(from = 1, to = framenumber), totaldisplacement[1:length(totaldisplacement[,4])-1,4])
  lines(seq(from = 1, to = framenumber), totaldisplacement[1:length(totaldisplacement[,4])-1,4])
  windows()
  plot(seq(from = 1, to = framenumber), totaldisplacement[1:length(totaldisplacement[,4])-1,5])
  lines(seq(from = 1, to = framenumber), totaldisplacement[1:length(totaldisplacement[,4])-1,5])
}
