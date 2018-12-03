#!/usr/bin/env Rscript

# Particle identification tracking and characterisation for 2D confocal images

confocaltrackroutine = function(remove_drift = TRUE){
  
  base_directory <- toString(getSrcDirectory(function(x) {x}))
  cat(base_directory, "\n")
  # source all scripts in the current directory
  setwd(base_directory)
  filelist <- c("pre_tracking/feature.r", "characterisation/gr2d.r", "tracking/iantrack.r", "tracking/driftremoval.R", "characterisation/isf.r", "pre_tracking/lowpass.r", "characterisation/msd.r", "pre_tracking/pretrack.r", "characterisation/shift.r", "characterisation/overcirc.r", "characterisation/psi6series.r", "characterisation/psi6loc.r", "characterisation/drawvoronoi.r")
  sapply(filelist,source,.GlobalEnv)
  library(EBImage)

  # File directory variables
  file_name <- "./documentation/brightfield"
  file_suffix = ".tif"
  num_digits = 4
  
  # put slash on end of dirname
  directory_name <- "./sample/"
  
  # Pretrack variables
  num_images <- 3             # How many images to read from file_name
  diameter <- 9               # Particle diameter - used in particle identification
  filter <- 9                 # Parameter for lowpass filter
  bgavg <- 9                  # Parameter for lowpass filter
  masscut <- 0.5              # Lowest integrated brightness for particle
  minimum <- 0.1              # Lowest pixel value for center of particle
  
  # Track variables
  edge_cutoff <- 20           # Cuts off this many pixels from each edge of the image in all data output - this is because particle identification is bad around the edges.
  max_displacement <- 10      # Used in tracking - the maximum allowed interframe displacement
  
  # Misc variables
  real_particle_size = 9      # Used as the particle size when doing characterisation
  timestep = 0.0294           # Frame time in seconds. Used for all data output to correct time in frames to time in seconds.
  gr_frames = 2               # How many frames of data to analyse for the g(r)
  trajectory_minimum = 1      # The smallest number of frames a trajectory must be linked for such that it is ouptut
  
  
  ### Main ###
  testimg <- readImage(paste(file_name, formatC(0, flag="0", width=num_digits), file_suffix, sep=""))
  
  image_x <- dim(testimg)[1]          #Width of the image in pixels
  image_y <- dim(testimg)[2]          #Height of the image in pixels
  
  # Identify particles in images
  pt <- pretrack(file_name, num_images, diameter, filter, bgavg, masscut, minimum, file_suffix, num_digits)
  ptfilt <- which(pt[,1] > edge_cutoff & pt[,1] < (image_x-edge_cutoff) & pt[,2] > edge_cutoff & pt[,2] < (image_y-edge_cutoff))
  
  write.table(pt, file=paste(directory_name, "raw_coords.txt", sep=""), row.names=FALSE, col.names=FALSE)

  #Print out the first image with tracked circles so the user can check the tracking is OK.  
  print_overcirc_image(file_name, file_suffix, pt, ptfilt, num_digits)

  calculate_num_particles(num_images, pt, ptfilt, timestep, directory_name)
  calcuate_gr(pt[ptfilt, ], edge_cutoff, image_x, image_y, gr_frames, real_particle_size, directory_name)
  calculate_psi_6(pt[ptfilt, ], directory_name, image_x, image_y, draw_voronoi=FALSE)
  
  tr <- do_tracking(pt, ptfilt, max_displacement, image_x, image_y, trajectory_minimum, remove_drift)
  
  write.table(tr, file=paste(directory_name, "tracked_coords.txt", sep=""), row.names=FALSE, col.names=FALSE)
  calculate_msd(tr, timestep, real_particle_size, directory_name) 
  calculate_isf(tr, real_particle_size, timestep, directory_name)
}


print_overcirc_image = function(file_name, file_suffix, pt, ptfilt, num_digits) {
  image_number = 0
  
  image_file_name <- paste(file_name,formatC(image_number, flag="0", width=num_digits), file_suffix , sep="") 
  imagecn <- channel(readImage(image_file_name), "gray")
  
  temp <- pt[ptfilt,]
  ptmask <- which(temp[, 6] == image_number)
  circledimage <- overcirc(imagecn, temp[ptmask,], rad=5) 
  #writeImage(circledimage, paste(directory_name, "overcirc", formatC(image_number, flag="0", width=num_digits), file_suffix, sep=""))
  
}

do_tracking = function(pt, ptfilt, max_displacement, image_x, image_y, trajectory_minimum, remove_drift) {
  # Track trajectories of identified particles
  if (max(pt[, 6]) < 2) {
    error("Cannot track particles if there are less than 3 frames")
  }
  tr <- iantrack(pretrack=pt[ptfilt,], maxdisp=max_displacement, imgsize=c(image_x,image_y), goodenough=trajectory_minimum)
  if (length(tr) == 0) {
    stop("No particles tracked")
  }
  
  if(remove_drift == TRUE) {
    tr <- driftremoval(tr, display_graphs = FALSE)
  }
  return(tr)
}

calculate_num_particles = function(num_images, pt, ptfilt, timestep, directory_name) {
  #Get particle count for each image
  particlecount <- matrix(nrow=num_images, ncol=1)
  particle_column <- pt[ptfilt,6]
  for(i in 1:num_images) {
    particlecount[i] <- sum(particle_column == i-1)
  }
  particlecount <- cbind(seq(1,num_images), seq(1,num_images)*timestep, matrix(particlecount))
  write.table(particlecount, file=paste(directory_name, "particlecount.txt", sep=""), row.names=FALSE, col.names=c("Frames", "Time", "Particle Count"))
}

calculate_msd = function(tr, timestep, particlesize, directory_name) {
  #Do msd and write it out
  if (max(tr[,6]) < 2) {
    stop("Cannot do msd for less than 3 frames")
  }
  msq <- msd(tr)
  msq <- cbind(msq[,1], msq[,1]*timestep, msq[,2:6], msq[,6]/particlesize, msq[,7])
  write.table(msq, file=paste(directory_name, "msd.txt", sep=""), row.names=FALSE, 
              col.names=c("Frames", "Time", "Mean dx", "Mean dy", "Mean dx^2", "Mean dy^2", "Mean dr^2 (Pixels)", "Mean dr^2 (Diameters)", "Particle Count"))  
}

calculate_isf = function(tr, particlesize, timestep, directory_name) {
  #Do isf and write it out
  if (max(tr[, 6]) < 2) {
    stop("Cannot do isf for less than 3 frames")
  }
  fsqt <- isf(tr, length=particlesize)
  fsqt <- cbind(fsqt[,1], fsqt[,1]*timestep, fsqt[,2:5])
  write.table(fsqt, file=paste(directory_name, "isf.txt", sep=""), row.names=FALSE, col.names=c("Frames", "Time", "Real", "Imaginary", "Modulus", "Samples"))
}

calcuate_gr = function (data, edge_cutoff, image_x, image_y, vargofrframes, particlesize, directory_name) {
  #Do g(r) and write it out
  if (vargofrframes > (max(data[, 6]) + 1)) {
    vargofrframes = (max(data[, 6]) + 1)
  }
  gr <- gr2d(data, nbins=200, deltar=0.5, imgsize=c(edge_cutoff, edge_cutoff, image_x-edge_cutoff, image_y-edge_cutoff), num_frames=vargofrframes)
  gr <- cbind(gr[,1], gr[,1]/particlesize, gr[,2]) # multiply the first column by timestep to give time in seconds
  write.table(gr, file=paste(directory_name, "gr.txt", sep=""), row.names=FALSE, col.names=c("Distance (Pixels)", "Distance (Diamaters)", "g(r)"))
}

calculate_psi_6 = function(tr, directory_name, image_x, image_y, draw_voronoi=FALSE) { 
  #Do psi_6
  psi6 <- psi6series(tr)
  write.table(psi6, file=paste(directory_name, "psi6.txt", sep=""), row.names=FALSE, 
              col.names=c("xcoords", "ycoords", "psi_6 real", "psi_6 imaginary", "psi_6 modulus", "frame number"))
  
  #Print the modulus
  cat("The mean modulus of psi6 is: ", mean(psi6[2:nrow(psi6),5]), "\n")

  if (draw_voronoi == TRUE) {
    vor <- drawvoronoi(psi6[2:nrow(psi6), 1], psi6[2:nrow(psi6), 2], psi6[2:nrow(psi6), 5], imgsize=c(0,image_x,0,image_y))
    writeImage(vor,file=paste(directory_name,"voronoi"), ".tif", sep="")
  }
}
