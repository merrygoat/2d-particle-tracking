# gr2d.r
# Ian Williams
# December 2013
#
# g(r) for two-dimensional data
# Taking into account the finite nature of microscope image
# Pretty slow right now

gr2d = function(data,nbins,deltar,imgsize,binary=FALSE,nframes=-1){
  
  # We are working with multiple frames - how many?
  # Note first frame is labeled zero, so add 1
  if (nframes==-1) {
    nframes <- max(data[,6]) + 1
  }
  
  # Output is 2 column matrix with r and g(r)
  # Or if binary is true we have 4 columns to hold bigbig, smallsmall and bigsmall
  if (binary==FALSE){
    output <- matrix(ncol=2,nrow=nbins)
    # Populate first column of output
    output[,1] <- seq(from=deltar,to=(nbins*deltar),by=deltar)
    output[,2] <- 0
  } else{
    output <- matrix(ncol=4,nrow=nbins)
    # Populate first column of output
    output[,1] <- seq(from=deltar,to=(nbins*deltar),by=deltar)
    output[,2] <- 0
    output[,3] <- 0
    output[,4] <- 0
  }
  
  
  # Need to sort density for binary case
  if(binary==FALSE){
    density <- nrow(data)/(((imgsize[3]-imgsize[1])*(imgsize[4]-imgsize[2]))*nframes)
  } else{
    w <- which(data[,6]==1)
    densitysmall <- length(w)/(((imgsize[3]-imgsize[1])*(imgsize[4]-imgsize[2]))*nframes)
    w <- which(data[,6]==0)
    densitybig <- length(w)/(((imgsize[3]-imgsize[1])*(imgsize[4]-imgsize[2]))*nframes)
    
    density <- nrow(data)/(((imgsize[3]-imgsize[1])*(imgsize[4]-imgsize[2]))*nframes)
  }
  
  #cat("Frames :",nframes,"\n")
  
  # Loops are nested thusly: frames, particles, bins
  # i labels frames,
  # j labels particles,
  # k labels bins
  
  cat("Progress of g(r)\n")
  objprogress <- txtProgressBar(min=0, max=nframes, style=3)
  
  for (i in 1:nframes){
    
    setTxtProgressBar(objprogress, i)  
    
    #cat("Frame ",i,"\n")
    
    # Get data for this frame
    wframe <- which(data[,6]==(i-1))
    thisframe <- data[wframe,]
    
    # How many particles in this frame?
    thisframeparticles <- nrow(thisframe)
    
    # Make array to hold this frame output
    if (binary==FALSE){
      thisframeout <- matrix(nrow=nbins,ncol=1)
      thisframeout[,1] <- 0
    } else{
      thisframeout <- matrix(nrow=nbins,ncol=3)
      thisframeout[,] <- 0
    }
    
    # Loop over particles
    for (j in 1:thisframeparticles){
      #cat("Particle :",j,"\n")
      # We are looking at particle j
      # How far are all the other particles from j?
      if (binary==FALSE){
        distfromj <- matrix(ncol=1,nrow=thisframeparticles)
        for (l in 1:thisframeparticles){
          distfromj[l,1] <- sqrt((thisframe[j,1] - thisframe[l,1])^2+(thisframe[j,2] - thisframe[l,2])^2)        
        }
      } else{
        # In binary case second column holds the bigbig smallsmall bigsmall identity
        distfromj <- matrix(ncol=2,nrow=thisframeparticles)
        for (l in 1:thisframeparticles){
          distfromj[l,1] <- sqrt((thisframe[j,1] - thisframe[l,1])^2+(thisframe[j,2] - thisframe[l,2])^2)        
          # 1-0 or 0-1 case which is bigsmall
          if (abs(thisframe[j,6] - thisframe[l,6])==1){distfromj[l,2] <- 2}
          # 1-1 case, which I think is smallsmall
          if ((thisframe[j,6]==1)&(thisframe[l,6]==1)){distfromj[l,2] <- 1}
          # 0-0 case, which I think is bigbig
          #cat("Here! \n")
          if ((thisframe[j,6]==0)&(thisframe[l,6]==0)){distfromj[l,2] <- 0}
        }
      }
      # distfromj holds the distance of all particles from j in single column
      
      # Loop over bins now
      for (k in 1:nbins){
        
        # Get start and end points for bin
        binend <- output[k,1]
        if (k==1){
          binstart <- 9E-9
        } else{
          binstart <- output[(k-1),1]
        }
        
        # Count how many particles are located between binstart and binend from particle j
        if (binary==FALSE){
          w <- which((distfromj[,1] > binstart)&(distfromj[,1] <= binend))
          if (is.numeric(length(w))){
            count <- length(w)
            #cat(count,"\n")
          } else{
            count <- 0
          } 
        } else{
          # Binary case
          # bigbig case
          w <- which((distfromj[,1] > binstart)&(distfromj[,1] <= binend)&(distfromj[,2]==0))
          if (is.numeric(length(w))){
            countbb <- length(w)
            #cat(countbb,"\n")
          } else{
            countbb <- 0
          } 
          # smallsmall case
          w <- which((distfromj[,1] > binstart)&(distfromj[,1] <= binend)&(distfromj[,2]==1))
          if (is.numeric(length(w))){
            countss <- length(w)
            #cat(countss,"\n")
          } else{
            countss <- 0
          } 
          # bigsmall
          w <- which((distfromj[,1] > binstart)&(distfromj[,1] <= binend)&(distfromj[,2]==2))
          if (is.numeric(length(w))){
            countbs <- length(w)
            #cat(countbs,"\n")
          } else{
            countbs <- 0
          } 
        }
      
        
        # Need to normalise this by the area of the annulus that falls within the image
        # First reset the edge-detection counters
        #edgexmin <- FALSE
        #edgexmax <- FALSE
        #edgeymin <- FALSE
        #edgeymax <- FALSE
        
        thetaxmin <- 0
        thetaymin <- 0
        thetaxmax <- 0
        thetaymax <- 0
                
        # Is particle j within distance binend of the edge of the image?
        if (thisframe[j,1] - imgsize[1] < binend){ 
          #edgexmin <- TRUE 
          thetaxmin <- 2*acos((thisframe[j,1] - imgsize[1]) / binend)
        }
        if (thisframe[j,2] - imgsize[2] < binend){
          #edgeymin <- TRUE 
          thetaymin <- 2*acos((thisframe[j,2] - imgsize[2]) / binend)
        }
        if ((imgsize[3] - thisframe[j,1]) < binend){ 
          #edgexmax <- TRUE 
          thetaxmax <- 2*acos((imgsize[3] - thisframe[j,1]) / binend)
        }
        if ((imgsize[4] - thisframe[j,2]) < binend){ 
          #edgeymax <- TRUE 
          thetaymax <- 2*acos((imgsize[4] - thisframe[j,2]) / binend)
        }
        
        # Now consider double counted angles for if more than 1 edge is crossed
        thetaxminymin <- (thetaxmin / 2) + (thetaymin / 2) - pi/2
        if (thetaxminymin <= 0){thetaxminymin <- 0}
        thetaxminymax <- (thetaxmin / 2) + (thetaymax / 2) - pi/2
        if (thetaxminymax <= 0){thetaxminymax <- 0}
        thetaxmaxymin <- (thetaxmax / 2) + (thetaymin / 2) - pi/2
        if (thetaxmaxymin <= 0){thetaxmaxymin <- 0}
        thetaxmaxymax <- (thetaxmax / 2) + (thetaymax / 2) - pi/2
        if (thetaxmaxymax <= 0){thetaxmaxymax <- 0}
        
        # Now find the angle of the arc at binend that actually sits in the image
        thetainimage <- (2*pi) - thetaxmin - thetaymin - thetaxmax - thetaymax + thetaxminymin + thetaxminymax + thetaxmaxymin + thetaxmaxymax
        
        # Normalise the count by the area
        if (binary==FALSE){
          count <- count / (binstart*thetainimage*deltar)
          thisframeout[k,1] <- thisframeout[k,1] + count
        } else{
          countbb <- countbb / (binstart*thetainimage*deltar)
          countss <- countss / (binstart*thetainimage*deltar)
          countbs <- countbs / (binstart*thetainimage*deltar)
          thisframeout[k,1] <- thisframeout[k,1] + countbb
          thisframeout[k,2] <- thisframeout[k,2] + countss
          thisframeout[k,3] <- thisframeout[k,3] + countbs
        }
      }
    }
    
    # Normalise by number of particles
    if (binary==FALSE){
      thisframeout[,1] <- thisframeout[,1] / (thisframeparticles*density)
      output[,2] <- output[,2] + thisframeout[,1]
    } else{
      thisframeout[,1] <- thisframeout[,1] / (thisframeparticles*densitybig)
      output[,2] <- output[,2] + thisframeout[,1]
      thisframeout[,2] <- thisframeout[,2] / (thisframeparticles*densitysmall)
      output[,3] <- output[,3] + thisframeout[,2]
      thisframeout[,3] <- thisframeout[,3] / (thisframeparticles*density)
      output[,4] <- output[,4] + thisframeout[,3]
    }
        
    plot(output[,1],output[,2],type="l")
    
  }
  
  close(objprogress)
  
  # Normalise by frames
  if (binary==FALSE){
    output[,2] <- output[,2] / nframes
  } else{
    output[,2] <- output[,2] / nframes
    output[,3] <- output[,3] / nframes
    output[,4] <- output[,4] / nframes
  }

  # Return output
  return(output)
  
}