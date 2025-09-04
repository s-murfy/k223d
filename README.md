
# k223d: generating stochastic slip distributions on unstructured mesh

<p align="center">
  <img src=".github/assets/SlipEg.png" alt="Overview of k223d" width="600"/>
</p>

[![DOI](https://zenodo.org/badge/111907406.svg)](https://doi.org/10.5281/zenodo.7525448)

## Overview

_k223d_ produces fractal stochastic slip distributions on non-planar faults that are described by triangular mesh. The programme will also provide the rupture time for each location on the fault plane for a given nucleation location and rupture velocity. 


In summary this progamme will:
- provide a stochastic slip distribution for each cell. 
- generate surface rupture if requested
- calculate the rupture time across the fault for a given nucleation location and rupture velocity. 

This programme is based on the [slipk2](https://github.com/andherit/slipk2) code for the generation of fractal slip while the kernal for calculating the distance across an unstructured mesh is based on [trilateration](https://github.com/andherit/trilateration)


➡️ For full documentation, see the [Wiki](https://github.com/s-murfy/k223d/wiki)


## Quick Start
Clone the repository:
```bash
git clone https://github.com/s-murfy/k223d.git
cd k223d
```




