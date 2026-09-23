# k223d Python examples

Examples here are organised by what they demonstrate, not by when they were added.

## `fault_sources/`
Tutorials showing how to get a fault mesh from a given geometry database and run k223d on it. These are meant to
be generic and reusable — the earthquake case studies below should link to these rather than re-explaining the
same load/remesh/run steps.

| Folder | Data source |
| --- | --- |
| `cfm_california/` | [SCEC Community Fault Model](https://www.scec.org/science/community-fault-model/) |
| `efsm_europe/` | [European Fault-Source Model 2020 (EFSM20)](https://seismofaults.eu/efsm20) |

Planned: a `slab2_subduction/` tutorial covering [Slab 2.0](https://www.usgs.gov/data/slab2-a-comprehensive-subduction-zone-geometry-model),
factored out of the mesh-building steps currently embedded in `earthquakes/kuril_2025/`.

## `earthquakes/`
Case studies of specific real earthquakes: taking a fault geometry (usually from one of the sources above),
conditioning slip on an observed/finite-fault model, and generating stochastic realisations with k223d. Each
folder's README should say which fault source it uses and link back to the relevant `fault_sources/` tutorial.

| Folder | Event | Fault source |
| --- | --- | --- |
| `kuril_2025/` | 2025 M8.8 Kuril earthquake | Slab 2.0 (mesh-building not yet split out, see its README) |

## Naming convention
- lowercase `snake_case`, no spaces
- `fault_sources/<source>_<region-or-scope>/` (e.g. `cfm_california`, `efsm_europe`, `slab2_subduction`)
- `earthquakes/<region>_<year>/` (e.g. `kuril_2025`)
- every example folder gets its own `README.md` describing the workflow and citing the data source(s) used

See also `../fortran/` for Fortran-only examples and `../local_py_scripts/` for shared helper functions imported
by these notebooks.
