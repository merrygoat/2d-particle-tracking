args <- commandArgs(trailingOnly=TRUE)
source("shift.r")
source("isf.r")

print(args)
l <- as.numeric(args[2])
cat("Loading data...\n")
data=read.table(args[1])
cat("First block...\n")
first <- isf(data,interval=1,length=l,maxtime=100)
cat("Second block...\n")
second <- isf(data,interval=100,length=l,maxtime=1000)
cat("Third block...\n")
third <- isf(data,interval=1000,length=l,maxtime=5000)
cat("Concatenating...\n")
alltogether <- rbind(first, second, third)
cat("Writing ISF...\n")
write.table(alltogether,paste("ISF",args[1], sep="_"))