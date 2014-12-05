# iantrack.r
# Feb 2014
#
# A new tracking routine that works differently
# Hopefully avoids the N^2 scaling in the old routine
#
# Takes pretrack data
# maxdisp is range to look for particle identifications
# imgsize is the x and y size of the original images in form c(x,y)
# goodenough is the minimum length of trajectories to return


iantrack = function(pretrack,maxdisp,imgsize,goodenough=0,memory=FALSE){
  
  
  nnewsum <- 0
  
  # We assume time starts at zero and increases to tmax
  tmax <- max(pretrack[,6])
  
  # Make output array of 6 columns to store:
  # x, y, mass, r^2, frame, particle number
  # 1 row for now, but we will ignore the first row at the end
  output <- matrix(nrow=1,ncol=7)
  #cat(output) #currently says NA NA NA NA NA NA
  
  # Label particles in first frame sequentially from 0
  # These form our initial particle identifications
  wzero <- which(pretrack[,6] == 0)
  zeropt <- pretrack[wzero,]
  zeroparticles <- nrow(zeropt)
  initident <- matrix(nrow=zeroparticles,ncol=1)
  initident[,1] <- seq(from=1,to=zeroparticles,by=1)
  output <- rbind(output,cbind(zeropt,initident))
  # So now particles in frame zero have an ID number 
  
  # Remove first row of output that contains a bunch of NA
  output <- output[2:(zeroparticles+1),]
  
  # The ID for the next new particle is...
  nextnewparticle <- zeroparticles + 1
  
  # Building circular mask of range maxdisp
  x2 <- seq(from=(0 - maxdisp),to=maxdisp,by=1)^2
  x2 <- t(x2)
  xmat <- matrix(nrow=((2*maxdisp) + 1),ncol=((2*maxdisp) + 1))
  xmat[,] <- x2
  ymat <- t(xmat)
  r2 <- xmat + ymat
  circmask <- r2<=maxdisp^2
  # circmask contains a logical matrix which is true within a circle of radius range
  # convert logical to numeric
  circmask <- circmask*1
  
  # I guess we'll have to loop over frames at some point
  # Start in second frame and compare it to the first
  for (i in 1:tmax){
    
    #cat("Tracking frame: ",i,"\n")
    
    # Get this frame co-ordinates
    wnow <- which(pretrack[,6] == i)
    nowcoords <- pretrack[wnow,]
    
    # How many particles in this frame?
    nowparticles <- nrow(nowcoords)
    
    # Get previous frame co-ordinates
    # Previous frame data comes from output array so it contains the particle identifications
    wbefore <- which(output[,6] == (i-1))
    beforecoords <- output[wbefore,]
    
    # How many particles in previous frame?
    beforeparticles <- nrow(beforecoords)
    
    # Make array nowident to hold particle IDs for nowparticles
    nowident <- matrix(nrow=nowparticles,ncol=1)
    # Make nowident -1 everywhere as a placeholder for no identification made
    nowident[,1] <- -1
    
    # Make "images" of the co-ords in both now and before
    # Put zero everywhere except for in co-ords nearest particle centres
    # In the particle locations put the row ID for the corresponding particle
    # We label rows first and cols second, so put x in rows and y in cols
    # This is technically transposed image
    # Add a border of size maxdisp around edges
    nowimage <- matrix(nrow=(imgsize[1]+(2*maxdisp)),ncol=(imgsize[2]+(2*maxdisp)))
    nowimage[,] <- 0
    for (j in 1:nowparticles){
      # Find nearest x and y pixel to particle j location
      nearx <- round(nowcoords[j,1]) + maxdisp
      neary <- round(nowcoords[j,2]) + maxdisp
      # We know the true co-ords are in row j of nowcoords
      # So label the nearest integer pixel in the image with the row number j
      nowimage[nearx,neary] <- j
    }
    # Actually I don't think I need beforimage
    #beforeimage <- matrix(nrow=imgsize[1],ncol=imgsize[2])
    #beforeimage[,] <- 0
    #for (j in 1:beforeparticles){
      # Find nearest x and y pixel to particle j location
      #nearx <- round(beforecoords[j,1])
      #neary <- round(beforecoords[j,2])
      # We know the true co-ords are in row j of nowcoords
      # So label the nearest integer pixel in the image with the row number j
      #beforeimage[nearx,neary] <- j
    #}
    
    # make xlo, xhi, ylo and yhi for beforecoords
    beforexlo <- beforecoords[,1]
    beforexhi <- beforecoords[,1] + (2*maxdisp)
    beforeylo <- beforecoords[,2]
    beforeyhi <- beforecoords[,2] + (2*maxdisp)
    # These are used to define the subimage window
    
    # Memory stuff
    # MEMORY DOESN'T WORK DON'T USE IT
    #if (i != 1){
     # wbeforebefore <- which(output[,5] == (i-2))
      #beforebeforecoords <- output[wbeforebefore,] 
      
      #beforebeforexlo <- beforebeforecoords[,1]
      #beforebeforexhi <- beforebeforecoords[,1] + (2*maxdisp)
      #beforebeforeylo <- beforebeforecoords[,2]
      #beforebeforeyhi <- beforebeforecoords[,2] + (2*maxdisp)
    #}
    
    # Loop over particles in the before image
    for (j in 1:beforeparticles){
      #cat("Looking at particle ",j,"\n")
      # Look in the now image in a subimage window located at the approx position of particle j in the before image
      nowsubimg <- nowimage[beforexlo[j]:beforexhi[j],beforeylo[j]:beforeyhi[j]]
      # Multiply this by circmask so we only look in a circular region
      
      #cat(nrow(nowsubimg)," ",ncol(nowsubimg),"\n")
      #cat(nrow(circmask)," ",ncol(circmask),"\n")
      nowsubimg <- nowsubimg*circmask
      
      # Where is this subimage non-zero (i.e. where are possible particles?)
      w <- which(nowsubimg != 0,arr.ind=TRUE)
      # How many possible particles?
      npossible <- nrow(w)
      #cat(npossible,"\n")
      # If there are no possible particles then this is non-numeric... This is a particle disappearance
      # If there is only one possible particle then this is our identification
      if (npossible==1){
        # Case of only one possibility
        # Particle must share ID with the only candidate
        id <- nowsubimg[w[1,1],w[1,2]]
        nowident[id,1] <- beforecoords[j,7]
        #cat(nowsubimg[w],"\n")
        #cat(id," ",beforecoords[j,6],"\n")
         
      } else if (npossible > 1){
        # How far is each possible particle from the centre of the subimage?
        #cat("Hello I am here \n")
        possibledisps <- matrix(nrow=npossible,ncol=2)
        for (k in 1:npossible){
          # This is the displacements loop
          # Location of particle j in before frame
          jx <- beforecoords[j,1]
          jy <- beforecoords[j,2]
          #cat(jx," ",jy,"\n")
                    
          # What particle number is possible k?
          kident <- nowsubimg[w[k,1],w[k,2]]
          #cat(j," ",kident,"\n")
          # Which has coords in now?
          kx <- nowcoords[kident,1]
          ky <- nowcoords[kident,2]
          #cat(kx," ",ky,"\n")
          
          # Find displacement between j in before and kident in now
          # first col of possibledisps holds the displacement
          possibledisps[k,1] <- sqrt((jx - kx)^2+(jy - ky)^2)
          # Second col holds the index ID of particle k
          possibledisps[k,2] <- kident
          #return(possibledisps)
                  
        }
        #cat(possibledisps,"\n")
        # Find the minimum displacement
        mindisp <- min(possibledisps[,1])
        # Which is it?
        wheremin <- which(possibledisps[,1]==mindisp)
        minid <- possibledisps[wheremin,2]
        # The row number in wheremin is then the identified particle
        nowident[minid,1] <- beforecoords[j,7]
      } else {
        # Particle disappearance handling here
        
        #if (memory == TRUE){
          # Now let's give the code a one-frame memory
          # Makes no sense in second frame
          #if (i != 1){
            #memsubimg <- nowimage[beforebeforexlo[j]:beforebeforexhi[j],beforebeforeylo[j]:beforebeforeyhi[j]]
            #memsubimg <- memsubimg*circmask
            
            #if (npossible==1){
              # Case of only one possibility
              # Particle must share ID with the only candidate
              #id <- memsubimg[w[1,1],w[1,2]]
              #nowident[id,1] <- beforebeforecoords[j,6]
              #cat(nowsubimg[w],"\n")
              #cat(id," ",beforecoords[j,6],"\n")
              
            #} else if (npossible > 1){
              # How far is each possible particle from the centre of the subimage?
              #cat("Hello I am here \n")
              #possibledisps <- matrix(nrow=npossible,ncol=2)
              #for (k in 1:npossible){
                # This is the displacements loop
                # Location of particle j in before frame
                #jx <- beforebeforecoords[j,1]
                #jy <- beforebeforecoords[j,2]
                #cat(jx," ",jy,"\n")
                
                # What particle number is possible k?
                #kident <- memsubimg[w[k,1],w[k,2]]
                #cat(j," ",kident,"\n")
                # Which has coords in now?
                #kx <- nowcoords[kident,1]
                #ky <- nowcoords[kident,2]
                #cat(kx," ",ky,"\n")
                
                # Find displacement between j in before and kident in now
                # first col of possibledisps holds the displacement
                #possibledisps[k,1] <- sqrt((jx - kx)^2+(jy - ky)^2)
                # Second col holds the index ID of particle k
                #possibledisps[k,2] <- kident
                #return(possibledisps)
                
              #}
              #cat(possibledisps,"\n")
              # Find the minimum displacement
              #mindisp <- min(possibledisps[,1])
              # Which is it?
             # wheremin <- which(possibledisps[,1]==mindisp)
              #minid <- possibledisps[wheremin,2]
              # The row number in wheremin is then the identified particle
              #nowident[minid,1] <- beforebeforecoords[j,6]
            #} else{
              # There is no else here unless I want to nest the memory further
              
            #}
            
          #}
        #}
        
        # When trajectory ends let's count how long it is and remove it if it is less than good enough
        # Which trajectory has ended?
        #if (goodenough != 0){
          #wended <- which(output[,6]==beforecoords[j,6])
         # nended <- length(wended)
          #if (nended < goodenough){
            #wkeep <- which(output[,6] != beforecoords[j,6])
            #output <- output[wkeep,]
            
            # Where is this subimage non-zero (i.e. where are possible particles?)
            #w <- which(nowsubimg != 0,arr.ind=TRUE)
            # How many possible particles?
            #npossible <- nrow(w)
            
          #}
        #}
      }
    }
    
    # Particles are identified
    # But what about new particles in frame now?
    # Loop over the elements of nowident
    nnew <- which(nowident[,1]==-1)
    #cat(length(nnew),"\n")
    nnewsum <- nnewsum + length(nnew)
    #cat("Cumulative New :", nnewsum,"\n")
    
    for (j in 1:nowparticles){
      if (nowident[j,1]==-1){
        # If particle has not been given a particle ID assign it the next available ID
        nowident[j,1] <- nextnewparticle
        nextnewparticle <- nextnewparticle + 1
      }
    }
    
    # Check that identifications are unique
    unq <- unique(nowident[,1])
    nunq <- length(unq)
    
    #cat(nunq," ",nrow(nowident),"\n")
    if (nunq==nowparticles){
      cat("Frame ",i,": Uniquely identified all particles.\n")
    } else{
      cat("Frame ",i,": Non-unique particles IDs. Problems!\n")
    }
    
    # Bind identifications to co-ords and then put this on the bottom of output
    nowout <- cbind(nowcoords,nowident)
    output <- rbind(output,nowout)
      
  }
  
  cat("Unique trajectories: ",length(unique(output[,7])),"\n")
  
  # If we have set the goodenough parameter let's get rid of trajectories that are too short
  if (goodenough != 0){
    len <- tabulate(output[,7])
    w <- which(len >= goodenough)
    ids <- seq(from=1, to=max(output[,7]), by=1)
    ids <- ids[w]
    output <- output[output[,7] %in% ids,]
    # Now relabel particles as a sequence starting from 1 to avoid gaps in particle labelling
    output[,7]<-as.numeric(factor(output[,7], levels = unique(output[,7])))
    cat("Sufficiently long trajectories: ",length(unique(output[,7])),"\n")
  }
  
  
  return(output)
  
}