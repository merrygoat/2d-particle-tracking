# chi4.r
# Peter Crowther
# computes chi 4 from tracked data
# July 2015
# Does not currently work with binary mixtures

chi4 = function(data){
  
  # We are working with multiple frames - how many?
  # Note first frame is labeled zero, so add 1
  nframes <- max(data[,6]) + 1
  
  # How many trajectories are we working with?
  # Note first trajectory is labeled one so don't add 1
  ntraj <- max(data[,7])
  
  # Loops are nested thusly: time, reference particle, comparing particle
  # reftime labels reference time (t0)
  # curtime labels current time (t),
  # j labels reference particle,
  # k labels comparing particle
  
  reftime <- 0    #for now we leave reftime = 0. Reference frame is first frame only.
  # frames are labeled from zero but start in 1 since 0 is the reference frame.

  # Get data for the ref frame 
  refframemask <- which(data[,6] == reftime)     # data mask
  refframe <- data[refframemask,]     # actual data
  refparticlelist <- c(refframe[,7])  # list of particle ID's in this frame
  
  for (curtime in reftime+1:nframes){
    cat("Time ",curtime,"\n")
    
    childframe <- c()
    # Get a list of particles in frame curtime that exist in refframe
    for (i in 1:length(refparticlelist)) {
      childframemask <- which(data[,6] == curtime & data[,7] == refparticlelist[i])
      childframe = rbind(childframe, data[childframemask,]) # iterate through and add each child particle to the list
    }
    
    browser()
    distfromj <- matrix(ncol=1,nrow=length(refparticlelist))
    for (j in 1:length(refparticlelist)) {
      distfromj[j] = 0
      for (k in j+1:length(childframe)) {
        distfromj[j] <- distfromj[j] + (sqrt((refframe[j,1] - childframe[k,1])^2+(refframe[j,2] - childframe[k,2])^2))    
      }
    }