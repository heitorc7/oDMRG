# oDMRG

# Introduction

 Code repository for the ITensor based Open DMRG algorithm as explained in [Analysis of a density matrix renormalization group approach for transport in open quantum system](https://arxiv.org/URL).

# Requirements

 This code is written in cpp to be used together with the ITensor library (version 2.x). You should have it installed, as well as the needed dependencies for ITensor. For more information go to the [ITensor page](https://www.itensor.org/index.html)

# Files

 The bib folder contains all the functions used in the oDMRG, including
* the LRN-sites.h which can be used to create a LRN-type MPS,
* Liouv.cpp, which can be used to create the Liovillian MPO object,
* MakeIvec.cpp which can be used to set a initial state to be the Identity, which as shown in the paper, speeds up the calculations,
* WarmUp.cpp function for the warm up (see paper Sec. 5.5)
* MagneticProfile.cpp to easily create different magnetic fields to apply to your systems.

 The file oDMRG.cc is the main file, and it can be altered or it can be run by inputting, in this order
1. number of sites N
1. coupling strength gamma
1. first bath parameter f1
1. last bath parameter fN
1. magnetic field h
1. anisotropy Delta


# Simple usage

## Running oDMRG

After going into the folder and editing the make file to fit your system type in make, then test the file by running

##### ./oDMRG 10 0.1 0.8 0.2 0 1

To simulate 10 sites with coupling 0.1, baths at initial 0.8 / final 0.2 with no magnetic field (0) for the XXX model (Delta = 1).
Additionally, run the included script to export the simulation output to a text file instead of having it shown on the terminal.

## Plotting

If you have Mathematica, use the included DataPlot.nb notebook to create scripts, properly process the data from the simulation, and also plot it. It should be easy to use it.
The other notebook, bibFunc.nb contains the functions you can use to process the data.
