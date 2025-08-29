
# k223d: generating stochastic slip distributions on unstructured mesh

[![DOI](https://zenodo.org/badge/111907406.svg)](https://doi.org/10.5281/zenodo.7525448)



## Description of k223d

_k223d_ produces fractal stochastic slip distributions on non-planar faults that are described by triangular mesh. The programme will also provide the rupture time for each location on the fault plane for a given nucleation location and rupture velocity. It is possible to assign surface nodes 

_k223d_ is based on the composite source model technique ([Zeng et al. BSSA, 1994](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/94GL00367)) whereby a set of hierarchical circular patches of slip are randomly placed on the fault plane based on a spatial probablity density function. Each individual patch of slip has a spatial distribution defined as a circular crack ([Ruiz et al., GJI, 2011](https://academic.oup.com/gji/article/186/1/226/697919)). At each point on the fault plane the contribution from the patches of slip at that location are summed together leading to the production of a final slip distribution that is self-similar. A further development has been the inclusion of surface rupture in _k223d_ when surface nodes are specified. In this newest version it is now possible to calculate the travel time of the rupture front across the fault plane using a user defined rupture velocity and a nucleation location.  

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

### Compile k223d
_k223d_ is written in fortran and there needs to be compiled with a fortran compiler. In the case of the docker example we will use gfortran whose libraries are preloaded in the docker environement. To compile the k223d, click on the terminal images in the notebook. And type the following commands:
```
cd codes 
make 
```
The _k223d_ should now compile without generating errors. 


### Running Notebooks
This involves two steps 
1. Compile _k223d_
2. Start a jupyter notebook

### Open notebook 
After compiling _k223d_ compiled the Jupyter notebooks can be opened. These are found in `notebook` folder. There are two examples provided: in the folder `planar_eg` the file _planar_mesh.ipynb_ is a workflow on how to produce a planar mesh, run _k223d_ and view the results. In the folder `EFSM_eg` the file _read_EFSM20.ipynb_ provides an example on how to read a mesh from the EFSM database, remesh it, run _k223d_, view the results and generate an file that can be viewed in QGIS. 



## References 
If using _k223d_ we would appreciate if you could reference the following article:
- Herrero, A. and Murphy, S., 2018. Self-similar slip distributions on irregular shaped faults. Geophysical Journal International, 213(3), pp.2060-2070.

In the case of surface rupture the following reference is relavent:
- Murphy, S. and Herrero, A., 2020. Surface rupture in stochastic slip models. Geophysical Journal International, 221(2), pp.1081-1089.

Should you end up using the notebooks in the production of publications, the following references may need to be considered: 
- _Gmsh_ has been used in the generation of meshes, check their website on how best to reference them(https://gmsh.info). For exampel the authors appreciate referencing the following paper: C. Geuzaine and J.-F. Remacle, Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering, Volume 79, Issue 11, pages 1309-1331, 2009.
- Graphics are produced using [Plotly](https://plotly.com/chart-studio-help/citations/)
- A number of other python libraries were used in certain notebooks such as [Meshio](https://zenodo.org/records/1288334), [GeoPandas](https://zenodo.org/records/3946761#.Xy24LC2ZPOQ), [PyProj4](https://zenodo.org/records/4571637).
- More fault meshes from the European fault-source model 2020 can be found here:  https://seismofaults.eu



## Version History 
**Version 2:** Inclusion of rupture time. Standardise input and output of the mesh to a vtk file format. 

**Version 1:** Provide static slip distribution. 


