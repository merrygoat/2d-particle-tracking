# drawvoronoi.R
# Ian Williams
# December 2013
#
# Draws voronoi polygons in a nicer way than plot.voronoi (hopefully)
# Colouring = 1 colours by neighbours
# Colouring = 2 colours by phase (require psi6)
# Colouring = 3 colours by polycrystal grain ID (require polyx)
# Colouring = 4 colours by magnitude of psi_6 (require psi6)
# Spots = 0 is no spots
# Spots = 1 is spots by phase
# Spots = 2 is spots by psi_6
# Arrows = 1 draws arrows
# Arrows = 0 does not draw arrows

drawvoronoi = function(xcoords,ycoords,psi6,polyx,imgsize,inclusions=0,colouring=1,arrows=0,spots=0,bigcut=NA){
  
  library(tripack)
  
  vor <- voronoi.mosaic(xcoords,ycoords)
  plot.voronoi(vor,xlim=c(imgsize[1],imgsize[2]),ylim=c(imgsize[3],imgsize[4]),do.points=FALSE,main="",sub="",asp=1)
  
  # And extract the polygons
  polygons <- voronoi.polygons(vor)
  # How many polygons do we have?
  npolygons <- length(polygons)
  
  # Generate rainbow colour vector
  rbow <- rainbow(100)
  
  #cat(colour,"\n")
  
 
  
  
  # We will loop over polygons and calculate the properties of each in turn
  for (i in 1:npolygons){
    
    #cat(i,"\n")
    
    # Get the ith polygon
    thispolygon <- polygons[i]
    thispolygon <- unlist(thispolygon)
    
    # This is a vector of first all x then all y
    # Want to make a 2 column matrix
    nvertices <- length(thispolygon)/2
    
    thispolygonmatrix <- matrix(ncol=2,nrow=nvertices)
    thispolygonmatrix[,1] <- thispolygon[1:nvertices]
    thispolygonmatrix[,2] <- thispolygon[(nvertices+1):(2*nvertices)]
    
    # Now we have x co-ords of vertices in 1st column and y co-ords in 2nd
    # Vertices are labeled in anticlockwise order, which is good for me
    
    # If any of the vertices fall outside of the image then reject the polygon
    reject <- FALSE
    if (any(thispolygonmatrix[,1] < imgsize[1])){reject <- TRUE}
    if (any(thispolygonmatrix[,1] > imgsize[2])){reject <- TRUE}
    if (any(thispolygonmatrix[,2] < imgsize[3])){reject <- TRUE}
    if (any(thispolygonmatrix[,2] > imgsize[4])){reject <- TRUE}
    
    # Proceed if we do not want to reject the polygon
    if (reject==FALSE){
      # Select colour
      if (colouring==1){
        
        if (nvertices==6){colour<-"khaki1"}
        if (nvertices==5){colour<-"indianred1"}
        if (nvertices==4){colour<-"red3"}
        if (nvertices==7){colour<-"cornflowerblue"}
        if (nvertices==8){colour<-"royalblue4"}
        if (nvertices==9){colour<-"green"}
        if (nvertices==10){colour<-"white"}
        if (nvertices==11){colour<-"black"}
        if (nvertices==12){colour<-"purple"}
        
        areasum <- 0
        xcsum <- 0
        ycsum <- 0
        
        for (j in 1:nvertices){
          # set k = j+1, except when j=nvertices
          if (j == nvertices){k=1}else{k=j+1}
          
          # Let's do area first
          areasum <- areasum + ((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2]))
          xcsum <- xcsum + ((thispolygonmatrix[j,1] + thispolygonmatrix[k,1])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          ycsum <- ycsum + ((thispolygonmatrix[j,2] + thispolygonmatrix[k,2])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          
        }
        
        area <- areasum/2.0
        xc <- xcsum/(6.0*area)
        yc <- ycsum/(6.0*area)
        
        # Black for inclusions
        if (length(inclusions) > 1){
          #cat("Hello! \n")
          if (any((sqrt((xc - inclusions[,1])^2+(yc - inclusions[,2])^2)) < 1)){colour<-"gray48"}
          
        }
        
        # If bigcut is specified then draw large polygons as grey
        if (is.numeric(bigcut)){
          #areasum=0
          # Loop over vertices
          #for (j in 1:nvertices){
            # set k = j+1, except when j=nvertices
            #if (j == nvertices){k=1}else{k=j+1}
          
            # Let's do area first
            #areasum <- areasum + ((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2]))
            #xcsum <- xcsum + ((thispolygonmatrix[j,1] + thispolygonmatrix[k,1])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
            #ycsum <- ycsum + ((thispolygonmatrix[j,2] + thispolygonmatrix[k,2])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
            
          #}
                  
          # Now find the area and centroids
          #area <- areasum/2.0
          #xc <- xcsum/(6.0*area)
          #yc <- ycsum/(6.0*area)
          
          # Grey for large particles
          if (area > bigcut){colour <- "#A9A9A9"}
          
          # Black for inclusions
          #if (inclusions != 0){
            
            #if (any((sqrt((xc - inclusions[,1])^2+(yc - inclusions[,2])^2)) < 10)){colour<-"black"}
            
          #}
          
        }
      
        polygon(x=thispolygonmatrix[,1],y=thispolygonmatrix[,2],col=colour)
        
      } else if(colouring==2){
        
        # phase colouring
        
        areasum <- 0
        xcsum <- 0
        ycsum <- 0
        
        for (j in 1:nvertices){
          # set k = j+1, except when j=nvertices
          if (j == nvertices){k=1}else{k=j+1}
          
          # Let's do area first
          areasum <- areasum + ((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2]))
          xcsum <- xcsum + ((thispolygonmatrix[j,1] + thispolygonmatrix[k,1])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          ycsum <- ycsum + ((thispolygonmatrix[j,2] + thispolygonmatrix[k,2])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          
        }
        
        area <- areasum/2.0
        xc <- xcsum/(6.0*area)
        yc <- ycsum/(6.0*area)
        
        # Identify polygon with a given value of psi6
        
        if(nvertices==6){
        
          howfar <- sqrt((xc - psi6[,1])^2+(yc - psi6[,2])^2)
          
          psi6ref <- which(howfar==min(howfar),arr.ind=TRUE)
          #if(length(psi6ref)>2){psi6ref <- psi6ref[1,]}
          ni <- psi6ref[1]
          
          phase <- atan2(psi6[ni,4],psi6[ni,3])
          
          phase2 <- (phase + pi)/(2*pi)
          phase2 <- round(100*phase2)
          if (phase2==0){phase2 <- 100}
          colour <- rbow[(phase2)]
        }        
        
        if (nvertices==5){colour<-"gray50"}
        if (nvertices==4){colour<-"gray50"}
        if (nvertices==7){colour<-"gray50"}
        if (nvertices==8){colour<-"gray50"}
        
        # Black for inclusions
        if (length(inclusions) > 1){
          #cat("Hello! \n")
          if (any((sqrt((xc - inclusions[,1])^2+(yc - inclusions[,2])^2)) < 1)){colour<-"gray48"}
          
        }
        #cat(colour,"\n")
        polygon(x=thispolygonmatrix[,1],y=thispolygonmatrix[,2],col=colour)
        
      } else if(colouring==3){
        
        # ================================ #
        # polycrystalline grains colouring #
        # ================================ #
        
        areasum <- 0
        xcsum <- 0
        ycsum <- 0
        
        # how many grains?
        ngrains <- length(unique(polyx[,13]))
        
        for (j in 1:nvertices){
          # set k = j+1, except when j=nvertices
          if (j == nvertices){k=1}else{k=j+1}
          
          # Let's do area first
          areasum <- areasum + ((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2]))
          xcsum <- xcsum + ((thispolygonmatrix[j,1] + thispolygonmatrix[k,1])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          ycsum <- ycsum + ((thispolygonmatrix[j,2] + thispolygonmatrix[k,2])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          
        }
        
        area <- areasum/2.0
        xc <- xcsum/(6.0*area)
        yc <- ycsum/(6.0*area)
        
        # get howfar from the polyx locations
        howfar <- sqrt((xc - polyx[,1])^2+(yc - polyx[,2])^2)
        
        # Only proceed if the nearest particle in polyx is nearer than some cut off (5 for now)
        if (min(howfar)<5){
          # Get the row in polyx that is nearest to the current polygon centre
          grainrefrow <- which(howfar==min(howfar),arr.ind=TRUE)
          ni <- grainrefrow[1]
          grainnum <- polyx[ni,13]
          if ((grainnum%%2)==0){
            graincol <- round((grainnum/ngrains)*100)
          } else{
            graincol <- round(100-((grainnum/ngrains)*100))
          }
          
          colour <- rbow[(graincol)]
        } else{
          colour <- "gray75"
        }
        
        polygon(x=thispolygonmatrix[,1],y=thispolygonmatrix[,2],col=colour)
        
      } else if (colouring==4){
        
        # psi_6 colouring
        
        rbow <- rainbow(100,start=0,end=(4/6))
        
        areasum <- 0
        xcsum <- 0
        ycsum <- 0
        
        for (j in 1:nvertices){
          # set k = j+1, except when j=nvertices
          if (j == nvertices){k=1}else{k=j+1}
          
          # Let's do area first
          areasum <- areasum + ((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2]))
          xcsum <- xcsum + ((thispolygonmatrix[j,1] + thispolygonmatrix[k,1])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          ycsum <- ycsum + ((thispolygonmatrix[j,2] + thispolygonmatrix[k,2])*((thispolygonmatrix[j,1]*thispolygonmatrix[k,2])-(thispolygonmatrix[k,1]*thispolygonmatrix[j,2])))
          
        }
        
        area <- areasum/2.0
        xc <- xcsum/(6.0*area)
        yc <- ycsum/(6.0*area)
        
        howfar <- sqrt((xc - psi6[,1])^2+(yc - psi6[,2])^2)
        
        psi6ref <- which(howfar==min(howfar),arr.ind=TRUE)
        #if(length(psi6ref)>2){psi6ref <- psi6ref[1,]}
        ni <- psi6ref[1]

        colour <- rbow[(psi6[ni,5])*100]
        
        polygon(x=thispolygonmatrix[,1],y=thispolygonmatrix[,2],col=colour)
        
      }
    }
    
  }
  
  # Draw spots if required
  if (spots != 0){
    nspots <- nrow(psi6)
    
    # spots = 1 means colour spots by phase of psi_6
    # spots = 2 means colour spots by psi_6
    if (spots == 1){
      for (i in 1:nspots){
        phase <- atan2(psi6[i,4],psi6[i,3])
        phase2 <- (phase + pi)/(2*pi)
        phase2 <- round(100*phase2)
        colour <- rbow[phase2]
        symbols(x=psi6[i,1],y=psi6[i,2],circles=6,bg=colour,add=TRUE,inches=FALSE)
      }
    } else if (spots == 2){
        for (i in 1:nspots){
          colour <- rbow[round(100-(50*psi6[i,5]))]
          symbols(x=psi6[i,1],y=psi6[i,2],circles=6,bg=colour,add=TRUE,inches=FALSE)
        }
    }
    
    # Draw additional spots on inclusions
    # Quick, dirty and ugly
    #doinclusions <- length(inclusions)
    #if (doinclusions > 1){
     # ninclusions <- nrow(inclusions)
      #for (i in 1:ninclusions){
        #symbols(x=inclusions[i,1],y=inclusions[i,2],circles = 4,bg="black",add=TRUE,inches=FALSE)
      #}
      
    #}
    
  }
  
  # Now draw arrows if specified
  if (arrows==1){
    narrows <- nrow(psi6)
  
    for (i in 1:narrows){
    
      # angle of arrow
      #phase <- atan2(psi6[i,4],psi6[i,3])
      #if (phase < 0){phase <- pi - abs(phase)}
      
      arrowlength=8
      
      phase <- atan2(psi6[i,4],psi6[i,3])
      phase2 <- (phase + pi)/(2*pi)
      phase2 <- round(100*phase2)
      colour <- rbow[(phase2+1)]
    
      xstart <- psi6[i,1] + ((arrowlength/2)*cos(phase))
      ystart <- psi6[i,2] + ((arrowlength/2)*sin(phase))
      xend <- psi6[i,1] - ((arrowlength/2)*cos(phase))
      yend <- psi6[i,2] - ((arrowlength/2)*sin(phase))
    
      arrows(x0=xstart,y0=ystart,x1=xend,y1=yend,lwd=2,col="black",length=0.05,code=2,angle=30)
      
      
    }
    
  }
  
  
}