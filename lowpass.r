# 2d spatial lowpass filter
# Ian Williams
# 21/09/2013

# lobject must be odd integer
# Also subtracts background and rescales brightnesses

lowpass = function(image,lobject,bgavg=5,hpass=0,background="mean"){
  
  # Lowpass filter
  lp <- makeBrush(lobject,shape="disc",step=FALSE)^2
  lp <- lp/sum(lp)
  
  l <- filter2(image,lp)
  
  #Highpass if desired
  if (hpass != 0){
    hp <-  matrix(1, nc=3, nr=3)
    hp[2,2] = hpass
    h = filter2(image,hp)
    
    g <- l-h
  
   } else{
    
    g <- l
    
  }
  
  # Local averaging for background
  bg <- matrix(1,nc=bgavg,nr=bgavg)
  bg[,] <- 1.0/(bgavg^2)
  b <- filter2(image,bg)
  
  g <- l - b
  
  # Subtract median i.e. background
  if (background=="mean"){
    g <- g - mean(g)
  } else{
    g <- g-median(g)
  }
      
  # Anything negative is set to zero
  black <- which(g<0,arr.ind=TRUE)
  g[black] <- 0
  
  # Now scale to full dynamic range 0 to 1
  scale <- 1/max(g)
  g <- g*scale
  
  return(g)
  
}