# k223d tutorial : how to read an EFSM20 fault and place a slip distribution on it

This is a tutorial on how to use [k223d](https://github.com/s-murfy/k223d) with a fault geometry taken from the
[European Fault-Source Model 2020 (EFSM20)](https://seismofaults.eu/efsm20). See notebook: `read_EFSM20.ipynb`

The notebook shows how to:
- import a fault mesh from the EFSM20 database (`ITCF02R.json`)
- remesh the fault plane at higher resolution in Gmsh
- set nucleation location, depth-dependent rupture velocity, and surface-rupture flags
- run k223d and inspect the slip distribution and rupture front
- write output for QGIS (geojson) and Paraview (vtk)

As a worked example, the tutorial reproduces a fault geometry associated with the M_w 6.5 Norcia earthquake
(30/10/2016, central Italy): https://terremoti.ingv.it/en/event/8863681

## Reference
EFSM20: https://seismofaults.eu/efsm20
