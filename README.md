


<div align="center">
<h1>WormAlgorithmBLBQ</h1>

 <img src="visualisations/torus_logo2.png" width="140" />

<p>Continuous-Time Quantum Monte Carlo Simulation of the Bilinear-Biquadratic Spin-1 Chain.</p>

<br />

[Visualization](#visualization) |
[Getting started](#getting-started)|
[Download Thesis](https://github.com/durnwalder) 
<br />
</div>
<hr />
This repository contains a Quantum Monte Carlo algorithm along with some utilities and a thesis describing the algorithm in detail. The algorithm simulates a generalization of the AKLT model, the bilinear-biquadratic spin-1 chain. The model is a striking example of the rich physics encountered in low-dimensional quantum spin systems. The various phases of the model were extensively studied with different analytical and numerical techniques. Most of the early numerical work on this model was done using exact diagonalization and, more recently, DMRG.  Our numerical study of the model relies on the worm algorithm in the path Integral representation, which is formulated in continous time. The worm algorithm was invented in the 1990s by Prokovief et al. It has since become one of the most popular methods to study bosonic lattice models and spin systems. 

# Visualization
<p align="left">
  <img src="visualisations/insert2.gif" width="400" />
</p>
The algorithm samples worldline configurations with three species of bosonic particles (red-green-blue) and two types of interactions.  This bosonic particle model is equivalent to a Spin-1 model. The two black squares are what we call worm head and tail, hence the name worm algorithm. The head performs a random walk in space time. When head and tail meet, the worm gets removed and reinserted at a random point. For a diagrammatic representation of the algorithm see and for a detailed discussion see thesis.pdf .    </figcaption>

# Getting started
WormAlgorithmBLBQ is based on the [ALPSCore](https://github.com/ALPSCore/ALPSCore)
library for physics simulations. Refer to [their website](http://alpscore.org/) for more detailed installation
instructions. ALPSCore has the following prerequisites:

  * a C++11-capable compiler.
  * CMake build system (version 3.1 or later)
  * Boost (version 1.56 or later)
  * HDF5 library 1.8.x (version 1.10 has a known problem)



## Building and installing ALPSCore
On HPC clusters, you may need to load the required modules first, e.g.  `module load gcc`.
1. Download: <br />
` git clone https://github.com/ALPSCore/ALPSCore.git` or <a href=" git clone https://github.com/ALPSCore/ALPSCore.git">download .zip file</a><br />

2. Build: <br />
`cd ALPSCore` <br />
`mkdir build` <br />
`cd build` <br />
`cmake .. -DALPS_INSTALL_EIGEN=yes -DCMAKE_INSTALL_PREFIX=/where/to/install/ALPSCore` <br />
`make `<br />
`make test`<br />

3. Install: <br />
`make install`<br /><br />  
## Building and running WormAlgorithmBLBQ


1. Download: <br />
`git clone https://github.com/durnwalder/WormAlgorithmBLBQ` or <a href="https://github.com/durnwalder/WormAlgorithmBLBQ/archive/master.zip">download .zip file</a><br />

1. Compile: <br />
`cd WormAlgorithmBLBQ` <br />
`cd src` <br />
`mkdir build` <br />
`cmake ..` <br />
`make`<br />

1. Run: <br />
`./qmc_worm ../parameter_files/aklt.ini`<br /><br />  

