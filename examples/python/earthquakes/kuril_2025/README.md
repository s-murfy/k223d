# k223d case study : 2025 M8.8 Kuril earthquake

This notebook generates stochastic slip distributions similar to the 2025 M8.8 Kuril earthquake on the
Kamchatka–Kuril subduction zone. See notebook: `kuril_mesh.ipynb`

The workflow:
- extracts the Kamchatka–Kuril interface from [Slab 2.0](https://www.usgs.gov/data/slab2-a-comprehensive-subduction-zone-geometry-model) and meshes it
- builds a slip probability density function (PDF) from the USGS
  [Finite Fault Slip Distribution](https://earthquake.usgs.gov/earthquakes/eventpage/us6000qw60/finite-fault) (v5)
- sets a depth-dependent rupture velocity and a nucleation location based on the actual hypocentre
- generates multiple stochastic slip distributions and rupture fronts with k223d, honouring the long-wavelength
  component of the USGS model while varying the short-wavelength (stochastic) contribution on each run
- converts output to geographic coordinates for [HySEA](https://edanya.uma.es/hysea/models/tsunami-hysea)
  tsunami simulation and to geojson for [QGIS](https://qgis.org/)

**Note:** this notebook currently bundles two things that are conceptually separate — (1) building a mesh from
Slab 2.0, which is generic to any subduction-zone event, and (2) the Kuril-specific PDF conditioning and slip
generation. Once a dedicated Slab 2.0 tutorial exists under `../../fault_sources/`, the mesh-building steps here
should be factored out and this notebook should reference that tutorial instead of repeating it.

## References
Goldberg, D. E., P. Koch, D. Melgar, S. Riquelme, and W. L. Yeck (2022). Beyond the teleseism: Introducing regional
seismic and geodetic data into routine USGS finite-fault modeling. Seismol. Res. Lett., 93(6), 3308–3323.
[doi](https://doi.org/10.1785/0220220047)

Hayes, G. P. (2018) Slab2 - A Comprehensive Subduction Zone Geometry Model: U.S. Geological Survey data release.
[doi](https://doi.org/10.5066/F7PV6JNV)

Hayes, G. P. (2017). The finite, kinematic rupture properties of great-sized earthquakes since 1990. Earth Planet.
Sci. Lett., 468, 94–100. [doi](https://doi.org/10.1016/j.epsl.2017.04.003)

Herrero, A. and Murphy, S. (2018). Self-similar slip distributions on irregular shaped faults. Geophysical Journal
International, 213(3), pp.2060-2070. [doi](https://doi.org/10.1093/gji/ggy104)

Murphy, S. and Herrero, A. (2020). Surface rupture in stochastic slip models. Geophysical Journal International,
221(2), pp.1081-1089. [doi](https://doi.org/10.1093/gji/ggaa055)
