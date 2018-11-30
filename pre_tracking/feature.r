# feature.r
# Identifies circular features
# Pretty much like Eric's and Paddy's 2d code
# Except obviously worse because I wrote it
# Ian Williams
# 22/09/2013

feature = function(image,diameter,separation=0,masscut,minimum,iterate=FALSE,ecccut=1,verbose=0){
                     
  # Require odd diameter
  if ((diameter %% 2)==0){
    cat('Requires odd diameter. Adding 1...\n')
    diameter <- diameter + 1    
  }
  
  # x and y sizes of image
  nx <- nrow(image)
  ny <- ncol(image)
  
  #cat(nx,'\n')
  #cat(ny,'\n')
  
  # Put a black border around the image 
  imgborder <- matrix(nrow=(nx + 2*diameter),ncol=(ny + 2*diameter))
  imgborder[,] <- 0
  imgborder[(diameter+1):(diameter+nx),(diameter+1):(diameter+ny)] <- image
  
  # Now we need to identify local maxima in brightness
  # Width of sampling region is diameter
  range <- (diameter-1)/2
  
  # Building circular mask
  x2 <- seq(from=(0 - range),to=range,by=1)^2
  x2 <- t(x2)
  xmat <- matrix(nrow=diameter,ncol=diameter)
  xmat[,] <- x2
  ymat <- t(xmat)
  r2 <- xmat + ymat
  circmask <- r2<=range^2
  # circmask contains a logical matrix which is true within a circle of radius range
  # convert logical to numeric
  circmask <- circmask*1
  
  # Theta mask
  x <- seq(from=(0 - range),to=range,by=1)
  y <- seq(from=(0 - range),to=range,by=1)
  theta <- matrix(nrow=diameter,ncol=diameter)
  for(j in 1:diameter){
    theta[j,] <- atan2(y[j],x)
  }
  
  # Sine and cosine masks
  sinmask <- sin(2*theta)*circmask
  cosmask <- cos(2*theta)*circmask
  
  # Dilate finds local maxima
  dilateimg <- dilate(imgborder,circmask)
  #dilateimg <- t(dilateimg)
  
  # Check for where dilated image = image and is greater than minimum
  # Returns logical matrix that is TRUE at local maxima i.e. particle centres
  # Recommend minimum = 0.5 to start with
  r <- ((dilateimg==imgborder)&(imgborder > minimum))
  r <- r*1
  
  # Now to get pixel co-ordinates out of this matrix image of particle centres
  coords <- which(r==1,arr.ind=TRUE)
    
  xcoords <- coords[,1]
  ycoords <- coords[,2]
  
  # Number of particles found
  nparticles <- length(coords[,1])
  if (verbose==1){
  cat('Initial particles: ',nparticles,'\n')
  }

  # Low and high x and y for each particle centre
  xlo <- xcoords-range
  xhi <- xcoords+range
  ylo <- ycoords-range
  yhi <- ycoords+range
  
  
  # New separation treatment Feb 2014 (prompted by Lizzie's work)
  cimgborder <- imgborder*r
  for (i in 1:nparticles){
    # Look at subregion of image
    b <- cimgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]]
    # And look in a circular region
    b <- b*circmask
    # Where in this region is the maximum?
    wmax <- which(b==max(b,na.rm=TRUE),arr.ind=TRUE)
    #cat(nrow(wmax),"\n")
    # Case where there is just one maximum pixel
    if (max(b,na.rm=TRUE) != 0){
      if (nrow(wmax) == 1){
        # Is this maximum at the centre?
        # If not then this pixel is not our guy
        #cat("Case 1","\n")
        #cat(wmax,"\n")
        if (((wmax[1,1] - (range+1))^2 + (wmax[1,2] - (range+1))^2) != 0){ cimgborder[xcoords[i],ycoords[i]] <- 0 }
      } else {
        # Case where there are multiple maxima
        # Are any at the centre?
        # If not then this pixel is not our guy
        #cat("Case 2","\n")
        if (all(((wmax[,1] - (range+1))^2 + (wmax[,2] - (range+1))^2) != 0)){ 
          cimgborder[xcoords[i],ycoords[i]] <- 0 
        } else{
          # But if one of them is the centre, how about we take the centre of mass of all the maxmima?
          avgmaxx <- round(mean(wmax[,1])) 
          avgmaxy <- round(mean(wmax[,2]))
          #cat(avgmaxx," ",avgmaxy,"\n")
          # Rounded to nearest pixel (in the subimage we are considering)
          # Set all the subimage in cimgborder to zero
          cimgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]] <- 0
          # Then set the average maximum position to 1
          cimgborder[(xlo[i] + avgmaxx),(ylo[i] + avgmaxy)] <- 1
        
        }
      }
    }
  }
    
  # What is left are valid particles now
  coords <- which(cimgborder != 0,arr.ind=TRUE)
  
  xcoords <- coords[,1]
  ycoords <- coords[,2]
  
  # Number of particles found
  nparticles <- length(coords[,1])
    if (verbose==1){
  cat('Particles after separation 1: ',nparticles,'\n')
  }

  # Low and high x and y for each particle centre
  xlo <- xcoords-range
  xhi <- xcoords+range
  ylo <- ycoords-range
  yhi <- ycoords+range
  
  
  # Now we are using our circular mask again
  # Estimating the brightness mass of each particle
  mass <- matrix(ncol=1,nrow=nparticles)
  for (i in 1:nparticles){
    # Sums brightness in circular region around each particle centre
    mass[i,1] <- sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*circmask)
  }
  
  # Now to apply the masscut
  # Retain only information for sufficiently bright particles
  masstrue <- which(mass >= masscut)
  xcoords <- xcoords[masstrue]
  ycoords <- ycoords[masstrue]
  mass <- mass[masstrue]
  xlo <- xlo[masstrue]
  xhi <- xhi[masstrue]
  ylo <- ylo[masstrue]
  yhi <- yhi[masstrue]
  nparticles <- length(xcoords)
  coords <- matrix(ncol=2,nrow=nparticles)
  coords[,1] <- xcoords
  coords[,2] <- ycoords

  # Number of particles remaining after masscut
  if (verbose==1){
  cat('Particles after masscut: ',nparticles,'\n')
  }

  # Separation is very slow for large numbers of particles 
  # Probably don't use this bit
  if (separation != 0){
    # Find separations between identified particles
    # Put these into nparticles by nparticles matrix
    partseps <- matrix(ncol=nparticles,nrow=nparticles)
    for (i in 1:nparticles){
      for (j in 1:nparticles){
        # Put distance between particles into partseps
        partseps[i,j] <- sqrt((coords[i,1]-coords[j,1])^2 + (coords[i,2]-coords[j,2])^2)
      }
    }
  
    # Make logical array for particles too close to one another
    # where too close is defined by separation
    tooclose <- which(partseps < separation,arr.ind=TRUE)
    # But remove those that correspond to diagonal elements
    # i.e. particle too close to self
    ntooclose <- length(tooclose[,1])
    # make keep logical matrix indicating which particles are kept and which are thrown away
    keep <- matrix(ncol=1,nrow=nparticles)
    # Initialise keep to true for all particles, update to false in forthcoming loop if particle is bad
    keep[,1] <- TRUE
    # Now loop over tooclose conflicts
    for (i in 1:ntooclose){
      if (tooclose[i,1] != tooclose[i,2]){
        # If two DIFFERENT particles are too close choose the brightest one as the true particle
        # Put equality in first case to avoid weird bug I found
        # This is an example of crappy programming - it throws away BOTH particles
        if (mass[tooclose[i,1]] >= mass[tooclose[i,2]]){keep[tooclose[i,2],1] <- FALSE}
        if (mass[tooclose[i,1]] < mass[tooclose[i,2]]){keep[tooclose[i,1],1] <- FALSE}
      }
    }
  
  
    # Throw away the particles that don't make the cut
    xcoords <- xcoords[keep]
    ycoords <- ycoords[keep]
    mass <- mass[keep]
    xlo <- xlo[keep]
    xhi <- xhi[keep]
    ylo <- ylo[keep]
    yhi <- yhi[keep]
    nparticles <- length(xcoords)
    coords <- matrix(ncol=2,nrow=nparticles)
    coords[,1] <- xcoords
    coords[,2] <- ycoords
  
    # Number of particles remaining after separation resolution
    if (verbose==1){
    cat('Particles after separation resolution: ',nparticles,'\n')
    }    
  }
    
  # Centroiding for true centres
  xcentroids <- matrix(ncol=1,nrow=nparticles)
  ycentroids <- matrix(ncol=1,nrow=nparticles)
  # Masks for x and y displacement of pixels from particle centre
  x <- seq(from=(0 - range),to=range,by=1)
  #x <- seq(from=1,to=diameter,by=1)
  xmask <- matrix(nrow=diameter,ncol=diameter)
  xmask[,] <- x
  ymask <- t(xmask)
  ymask <- ymask*circmask
  xmask <- xmask*circmask
  
  for (i in 1:nparticles){
    # Sums brightness in circular region around each particle centre
    # x and y reverse again because we label rows first in R for some reason
    xcentroids[i,1] <- sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*xmask)
    ycentroids[i,1] <- sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*ymask)
  }
  
  xcentroids[,1] <- (xcentroids[,1]/mass[])
  ycentroids[,1] <- (ycentroids[,1]/mass[])
  
  # True centres here
  coords[,1] <- coords[,1] + xcentroids[,1]
  coords[,2] <- coords[,2] + ycentroids[,1]
  xcoords <- coords[,1]
  ycoords <- coords[,2]
  
  #cat(sum(xcentroids[,1])+sum(ycentroids[,1]),"\n")
  
  # Iteration of centroiding
  # Iteration counter
  # THIS STILL NEEDS SOME WORK
  itnum <- 0
  if (iterate==TRUE){
    repeat{
      # break if centroids are perfect
      if(abs(sum(xcentroids[,1])+sum(ycentroids[,1]))<1){break}
      # also break if we repeat too many times (for now let's say 10)
      if(itnum >= 20){break}
      
      # increment iteration number
      itnum <- itnum + 1
      
      xlo <- xcoords-range
      xhi <- xcoords+range
      ylo <- ycoords-range
      yhi <- ycoords+range
      
      for (i in 1:nparticles){
        # Sums brightness in circular region around each particle centre
        mass[i] <- sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*circmask)
      }
      
      for (i in 1:nparticles){
        # Sums brightness in circular region around each particle centre
        # x and y reverse again because we label rows first in R for some reason
        xcentroids[i,1] <- sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*xmask)
        ycentroids[i,1] <- sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*ymask)
      }
      
      xcentroids[,1] <- (xcentroids[,1]/mass[])
      ycentroids[,1] <- (ycentroids[,1]/mass[])
      
      # True centres here
      coords[,1] <- coords[,1] + xcentroids[,1]
      coords[,2] <- coords[,2] + ycentroids[,1]
      xcoords <- coords[,1]
      ycoords <- coords[,2]   
      
      cat(sum(xcentroids[,1])+sum(ycentroids[,1]),"\n")
      
    }
    
  }
  
  
  # Since we added a border of width diameter to the image, we should subtract this from the co-ords
  
  # Put stuff into output array
  output <- matrix(nrow=nparticles,ncol=5)
  output[,1] <- coords[,1] - diameter
  output[,2] <- coords[,2] - diameter
  output[,3] <- mass[]
  
  # r^2 mask for radius of gyration
  r2mask <- r2*circmask
  # Estimating radius of gyration and eccentricity
  rg <- matrix(ncol=1,nrow=nparticles)
  ecc <- matrix(ncol=1,nrow=nparticles)
  for (i in 1:nparticles){
    rg[i,1] <- sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*r2mask)
    ecc[i,1] <- sqrt(sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*cosmask)^2 + sum((imgborder[xlo[i]:xhi[i],ylo[i]:yhi[i]])*sinmask)^2)
  }
  rg[,1] <- rg[,1]/mass[]
  output[,4] <- rg[,1]
  ecc[,1] <- ecc[,1]/mass[]
  output[,5] <- ecc[,1]
  

  
  
  # Return the output array with data in columns:
  # x, y, mass, r_g, eccentricity

  w <- which(output[,5]<=ecccut)
  output <- output[w,]
   
  if (verbose==1){
  cat('Particles after eccentricity cut: ',nrow(output),'\n')
  }
  
  #Area fraction
  eta <- (pi/4 * diameter^2) * (nparticles/(ny*nx))
  cat('Area Fraction: ',eta,'\n')
  
  return(output)
  
  
  
}