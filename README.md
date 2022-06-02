


# WormAlgorithmBLBQ: Continuous-time quantum Monte Carlo simulation of the bilinear-biquadratic Spin-1 chain
| [Getting started](#getting-started)
| [The Worm Algorithm](#the-worm-algorithm)
| [Order Parameters](#order-parameters)

<p align="left">
  <img src="visualisations/insert.gif" width="400" />
</p>



- [WormAlgorithmBLBQ: Continuous-time quantum Monte Carlo simulation of the bilinear-biquadratic Spin-1 chain](#wormalgorithmblbq-continuous-time-quantum-monte-carlo-simulation-of-the-bilinear-biquadratic-spin-1-chain)
- [Getting started](#getting-started)
  - [Building and installing ALPSCore](#building-and-installing-alpscore)
  - [Building and running WormAlgorithmBLBQ](#building-and-running-wormalgorithmblbq)
  - [Parameter files](#parameter-files)
  - [Output files](#output-files)
- [The Worm Algorithm](#the-worm-algorithm)
  - [Negative-sign problem and model selection](#negative-sign-problem-and-model-selection)
  - [Requirements](#requirements)
  - [Worm type updates](#worm-type-updates)
- [Order Parameters](#order-parameters)
  - [The Haldane phase - String order parameter](#the-haldane-phase---string-order-parameter)
  - [Spin-Spin correlation function](#spin-spin-correlation-function)
  - [Dimer order parameter](#dimer-order-parameter)
  - [Off-diagonal long range order](#off-diagonal-long-range-order)





# Getting started
These codes are based on the [ALPSCore](https://github.com/ALPSCore/ALPSCore)
library. Refer to [their website](http://alpscore.org/) for installation
instructions. At the time of this writing, ALPSCore imposed the following system
requirements:

  * C++ compiler
  * CMake build system (version 3.1 or later)
  * Boost (version 1.56 or later)
  * HDF5 library 1.8.x (version 1.10 has a known problem)
  
Beyond these, our codes require

  * a C++11-capable compiler. We have only tested our code when ALPSCore is installed with modern compilers that have C++11 features turned on by default.


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
`git clone https://github.com/durnwalder/` or <a href="https://github.com/khavernathy/mcmd/archive/master.zip">download .zip file</a><br />

2. Compile: <br />
`cd mcmd` <br />
`cd src` <br />
`bash compile.sh   [ options ]` <br />
`cd ..` <br />
`export PATH=$PATH:/path/to/mcmd/`<br />

3. Run: <br />
`mcmd mcmd.inp`<br /><br />  


## Parameter files
<p align="center">
  <img src="visualisations/torus_logo2.png" width="150" />
</p>

## Output files

# The Worm Algorithm



The bilinear-biquadratic spin-1 chain is a striking example of the rich physics encountered in low-dimensional quantum spin systems. The various phases of the model were extensively studied with different analytical and numerical techniques. Most of the early numerical work on this model was done using exact diagonalization and, more recently, DMRG.  Our numerical study of the model relies on Quantum Monte Carlo, more specifically, the worm Algorithm in the path Integral representation. The worm algorithm was invented in the 1990s by Prokovief et al \cite{prokofev_exact_1998}. It is one of the most popular methods to study bosonic lattice models and spin systems. 

$$H=\sum_{i=1}^{N-1} \cos (\theta)  \vec{S}_{i} \otimes \vec{S}_{i+1}+\sin (\theta) \left(\vec{S}_{i} \otimes \vec{S}_{i+1}\right)^{2}+\sum_{i=1}^{N} D\left(S_{i}^{z}\right)^{2}$$

 ## Negative-sign problem and model selection   

Requirements
------------
 ## Worm type updates   

# Order Parameters
## The Haldane phase - String order parameter
<p>
  <img src="visualisations/worm_aklt.gif" width="400" />
</p>

## Spin-Spin correlation function
## Dimer order parameter
## Off-diagonal long range order


<p>
  <img src="visualisations/flow.png" width="500"/>
</p>
<p>
  <img src="visualisations/string.png" width="500"/>
</p>
