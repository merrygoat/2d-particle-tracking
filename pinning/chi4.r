# chi4.r
# Ian Williams
# May 2014

# Calculate (unormalised) chi_4 for 2d data as is my wont
# Cutoff is a length, typically 0.3 particle diameters (in pixels)

chi4 = function(tracks,cutoff){
  
  # Max frame number
  nframes <- max(tracks[,6])
  
  # Sort the tracks to be arranged by particle ID rather than frame
  tracks <- tracks[sort.list(tracks[,7]),]
  
  # Make log-spaced time intervals
  intervals <- round(1.15^(seq(from=1,to=100,by=1)))
  intervals <- unique(intervals)
  w <- which(intervals <= nframes)
  intervals <- intervals[w]
  
  # Make chi4 array
  chi4 <- matrix(nrow=length(intervals),ncol=2)
  chi4[,1] <- intervals
  
  # First we calculate Q(t), where t is the interval
  # Need to loop over intervals, then frames within the interval, I think
  for (i in 1:length(intervals)){
    # i labels interval
    
    cat("Interval :",i,"\n")
    
    qt <- qoft(tracks,dt=intervals[i],cutoff=cutoff)
    qt2 <- qt[,2]^2
    
    # Normalise chi4 properly in future    
    chi4[i,2] <- mean(qt2) - (mean(qt[,2]))^2
    
  }
  
  return(chi4)
  
  
}