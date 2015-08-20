# READING THINGS BACK IN
#####################################
setwd("Volumes/Seagate Backup Plus Drive/Projects/Pinning/Real_Experiments_Oct2014/JamesPeter/Peter2/Slide 5/2-55pm")

pt <- read.table("pretrack.dat",sep="\t")
pt <- data.matrix(pt)

tr <- read.table("track.dat",sep="\t")
tr <- data.matrix(tr)

trsize <- read.table("track_sizetag.dat",sep="\t")
trsize <- data.matrix(trsize)

pntag <- read.table("track_sizepintag.dat",sep="\t")
pntag <- data.matrix(pntag)

trbignonpin <- read.table("track_unpinned_big.dat",sep="\t")
trbignonpin <- data.matrix(trbignonpin)

trsmallnonpin <- read.table("track_unpinned_small.dat",sep="\t")
trsmallnonpin <- data.matrix(trsmallnonpin)

voro <- read.table("tr_voronoi.dat",sep="\t")
voro <- data.matrix(voro)

#####################################


# TESTING STUFF BEFORE BATCH TRACKING
#####################################
img <- channel(readImage("Images/sample5_2.55pm_00000001.tif"),"grey")
lp<-lowpass(img,lobject=11,bgavg=13)
f<-feature(lp,diameter=15,masscut=8,minimum=0.2)
#####################################

# ACTUAL TRACKING
#####################################
# binaryfeature.r is retired after cocking up one too many times
# Now the size identification is last thing we do after tracking!
#####################################

pt <- pretrack(filename="Images/sample5_2.55pm_0000",images=3740,startimg=1,diameter=15,masscut=8,minimum=0.2,filter=11,bgavg=13)
pt[,6] <- pt[,6] - 1
write(t(pt),file="pretrack.dat",ncolumns=6,sep="\t")

tr <- iantrack(pretrack=pt,maxdisp=5,imgsize=c(768,768),goodenough=2)
write(t(tr),file="track.dat",ncolumns=7,sep="\t")

# Cutoff between big and small particles. Do a histo here to determine the cutoff.
trsize <- sizetag(tr,cutoff=26)
write(t(trsize),file="track_sizetag.dat",ncolumns=8,sep="\t")
# 1 = big, 0 = small

# This is the old pintag method
#pntag <- pintagoct(trsize,longinterval=50,longthresh=1.0,posthresh=0.7)
# Remove the first 50 frames and renumber from 0
#w<-which(pntag[,6]>50)
#pntag <- pntag[w,]
#pntag[,6] <- pntag[,6]-50
#write(t(pntag),file="track_sizepintag.dat",ncolumns=9,sep="\t")

# New pintag based on sd position in time chunks
# Remember to alter the chunkwidth based on density
pntag <- sdpchpintag(trsize,chunkwidth=1000,sdthresh=0.3)

pntag <- pintidy(pntag)
write(t(pntag),file="track_sizepintag.dat",ncolumns=9,sep="\t")

pntime <- pinsintime(pntag)
write(t(pntime),file="pinsintime.dat",ncolumns=7,sep="\t")

particles <- binaryparticlecount(trsize)
write(t(particles),file="particlecount.dat",ncolumns=4,sep="\t")

# Split tracks up a bit
#######################

w <- which(pntag[,8]==1)
trbig <- pntag[w,]
write(t(trbig),file="track_big.dat",ncolumns=8,sep="\t")

w <- which(pntag[,8]==0)
trsmall <- pntag[w,]
write(t(trsmall),file="track_small.dat",ncolumns=8,sep="\t")

w <- which(pntag[,9]==1)
trpins <- pntag[w,]
write(t(trpins),file="track_pins.dat",ncolumns=9,sep="\t")

w <- which(pntag[,9]==0)
trnonpin <- pntag[w,]
write(t(trnonpin),file="track_unpinned.dat",ncolumns=9,sep="\t")

w <- which(trbig[,9]==0)
trbignonpin <- trbig[w,]
write(t(trbignonpin),file="track_unpinned_big.dat",ncolumns=9,sep="\t")

w <- which(trsmall[,9]==0)
trsmallnonpin <- trsmall[w,]
write(t(trsmallnonpin),file="track_unpinned_small.dat",ncolumns=9,sep="\t")




# Some static and dynamic analysis
###################################
w <- which(trsize[,6]<5)
gr <- gr2d(trsize[w,],nbins=500,deltar=0.1,imgsize=c(768,768),binary=TRUE)
write(t(gr),file="grbinary.dat",ncolumns=4,sep="\t")

msdbig <- msd(trbignonpin,drift=TRUE)
write(t(msdbig),file="msd_unpinned_big.dat",ncolumns=7,sep="\t")

msdsmall <- msd(trsmallnonpin,drift=TRUE)
write(t(msdsmall),file="msd_unpinned_small.dat",ncolumns=7,sep="\t")

isfbig <- isf(trbig,length=24.3,drift=TRUE)
write(t(isfbig),file="isf_all_big.dat",ncolumns=5,sep="\t")

isfsmall <- isf(trsmall,length=14.9,drift=TRUE)
write(t(isfsmall),file="isf_all_small.dat",ncolumns=5,sep="\t")


voro <- newvoroseries(trsize,imgsize=c(768,768),edgerm=50)
write(t(voro),file="tr_voronoi.dat",ncolumns=10,sep="\t")

vorphi <- voronoibinaryareafrac(voro)
write(t(vorphi),file="areafrac_voronoi.dat",ncolumns=8,sep="\t")

neigh <- neighboursbinary(voro)
write(t(neigh),file="binaryneighbourscount.dat",ncolumns=15,sep="\t")





