

k223d

Create a stochastic slip distribution with k2 amplitude spectra on an unstructured mesh

This program generates a k2 slip distribution on fault surface that is described by an unstructured (triangular mesh). A composite source model technique (Zeng et al. BSSA, 1994) is used to generated the self-similar slip distribution and the programme is based on the slipk2 code ( https://github.com/andherit/slipk2 ). The computation of the distance across the mesh used a double lateration technique (Herrero et al, in submission) the kernel of which can be found at https://github.com/andherit/trilateration )

The programme can be compiled using the 'install.sh' script, otherwise the makefile in the src directory can be used. The makefile and install.sh are currently set for a gfortran compiler. When compiled a k223d.x file is generated.

Usage: ./k223d.x in same folder as input_file and mesh file. The output is generated in the same folder as the programme is run and can be viewed using GMT or paraview (i.e. https://www.paraview.org ).


The mesh should contain only triangles (i.e unstructured) and be in abaqus format  (i.e. .inp ending ) - two examples are provided Scotia.inp and test.inp. To generate a planar fault with the program provided in the make_mesh folder - the output from this is 'test.inp' and can be used directly in k223d.  
