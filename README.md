# README - 2D Particle tracking code #

This code was developed in R by Ian Williams in 2014. It is a development of the code written by John Crocker and David Grier in IDL which is now hosted on Eric Weeks' website here: http://www.physics.emory.edu/faculty/weeks//idl/

The repository contains a PDF tutorial written by Ian that covers the basics of R as well as the specifics of compiling and running the code.

The code was substantially revised by Peter Crowther from 2015 to 2018.

## Installation 
The installation instructions in the Trackpack.pdf have been superceeded due to an updated version of the library. To install the EBImage and tripack libraries required to run the code use the commands:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage", version = "3.8")

install.packages("tripack")

## Running the code
Use the main.r script to do tracking. This can track a single or multiple frames.

Note that the overcirc function and g(r) have not been optimised for large numbers of particles > 1000. They will be very slow for larger numbers. For datasets with > 1000 particles it is a good idea to comment these lines out if they are not vital to your analysis.