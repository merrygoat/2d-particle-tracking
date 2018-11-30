# 2d spatial lowpass filter
# Ian Williams
# 21/09/2013

# lobject must be odd integer
# Also subtracts background and rescales brightnesses

lowpass = function(image,lobject,bgavg=5,background="mean"){
  
  # Subtract median i.e. background
  if (background=="mean"){
    image <- image - mean(image)
  } else if(background=="median"){
    image <- image-median(image)
  }
  
  # Anything negative is set to zero
  image[image < 0] <- 0
  
  # Now scale to full dynamic range 0 to 1
  image <- image/max(image)

  # Lowpass filter
  browser()
  lp_filter <- makeBrush(lobject,shape="disc",step=FALSE)^2
  lp_filter <- lp_filter/sum(lp_filter)
  
  # Filters the image using the fast 2D FFT convolution product
  l <- filter2(image, lp_filter)
  
  # Local averaging for background
  bg <- matrix(1, nc=bgavg, nr=bgavg)
  bg[,] <- 1.0/(bgavg^2)
  b <- filter2(image,bg)
  
  g <- l - b
  
  # Now image has been filtered, renormalise it
  g[g < 0] <- 0
  g <- g/(max(g))

  return(g)
  
}