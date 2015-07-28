#!/usr/bin/env Rscript

# source all scripts in the current directory
setwd("/Users/pc9836/Documents/git/2d-particle-tracking")
filelist <- c("chi4.r", "feature.r", "gr2d.r", "iantrack.r", "isf.r", "lowpass.r", "msd.r", "pretrack.r", "shift.r")
sapply(filelist,source,.GlobalEnv)
biocLite()
library(EBImage)

#Pretrack variables
varimages <- 200
vardiameter <-19
varfilter <-11
varbgavg <- 11
varmasscut <- 1
varminimum <- 0.1

#Track variables
vartrlengthscale <- 19
varimgx <- 400
varimgy <- 300
varedgecutoff <- 15
varmaxdisp <- 8

#More variables
varfilename <- "/Volumes/WIN_DATA/Confocal/STED/15-07-27/0.52-vi/FITC-0.52-vi-"
vardirname <- "/Volumes/WIN_DATA/Confocal/STED/15-07-27/0.52-vi/"

### Main ###

pt <- pt <- pretrack(filename=varfilename,images=varimages,diameter=vardiameter,filter=varfilter,bgavg=varbgavg, masscut=varmasscut,minimum=varminimum,chan="grey")
ptfilt <- which(pt[,1] > varedgecutoff & pt[,1] < varimgx-varedgecutoff & pt[,2] > varedgecutoff & pt[,2] < varimgy-varedgecutoff)

#Get particle count for each image
particlecount <- rep(0, varimages)
for(i in 0:varimages-1) {particlecount[i] <- sum(pt[ptfilt,6] == i)}
write(t(particlecount),file=paste(vardirname, "particlecount.txt"),ncolumns=1,sep="\t")

tr <- iantrack(pretrack=pt[ptfilt,],maxdisp=varmaxdisp,imgsize=c(varimgx,varimgy),goodenough=10)

#Don't need to filter tr as it was already filtered from pt
#Do msd and write it out
msq <- msd(tr)
write(t(msq),file=paste(vardirname, "msd.txt"),ncolumns=7,sep="\t")

#Do isf and write it out
fsqt <- isf(tr,length=19)
write(t(fsqt),file=paste(vardirname, "isf.txt"),ncolumns=5,sep="\t")

#Do g(r) and write it out
gr <- gr2d(tr[1:30000,],nbins=200,deltar=0.5,imgsize=c(varedgecutoff,varedgecutoff,varimgx-varedgecutoff,varimgy-varedgecutoff))
write(t(gr),file=paste(vardirname, "gr.txt"),ncolumns=2,sep="\t")