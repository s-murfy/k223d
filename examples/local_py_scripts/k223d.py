"""
k223d Python interface
======================
Wraps the Fortran k223d_compute_source C API via ctypes.

Requires libk223d.so in the same directory as this file.

Usage
-----
Slip only (no velocity or time in mesh file):

    import numpy as np
    import k223d

    nodes = np.loadtxt('nodes.txt')        # shape (nnodes, 3)
    cells = np.loadtxt('cells.txt', dtype=int)  # shape (ncells, 3), 1-based

    slip = k223d.compute_source(nodes, cells, mw=6.5)

Slip + rupture time:

    t0 = np.full(nnodes, 1.e32)
    t0[nucleation_node - 1] = 0.0          # nucleation node, 1-based index

    slip, time = k223d.compute_source(
        nodes, cells, mw=6.5,
        velocity=vel,
        time=t0,
    )

VTK field names expected by the Fortran reader: 'velocity', 'time', 'surface',
'pdf'. See README for mesh file format details.
"""

import numpy as np
import ctypes
import os
import warnings

# ---------------------------------------------------------------------------
# Load shared library
# ---------------------------------------------------------------------------
_lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libk223d.so')
try:
    _lib = ctypes.CDLL(_lib_path)
except OSError as e:
    raise ImportError(
        f"Could not load k223d shared library at {_lib_path}.\n"
        f"Build it with: make libk223d.so\n"
        f"Original error: {e}"
    )

# ---------------------------------------------------------------------------
# Ctypes pointer aliases
# ---------------------------------------------------------------------------
_dbl_p = ctypes.POINTER(ctypes.c_double)
_int_p = ctypes.POINTER(ctypes.c_int)

# ---------------------------------------------------------------------------
# Declare Fortran C API signature
# ---------------------------------------------------------------------------
_lib.k223d_compute_source.argtypes = [
    ctypes.c_int,    # nnodes
    ctypes.c_int,    # ncells
    _dbl_p,          # px[nnodes]
    _dbl_p,          # py[nnodes]
    _dbl_p,          # pz[nnodes]
    _int_p,          # cells[ncells,3] C-order == Fortran cells(3,ncells)
    ctypes.c_double, # mw
    ctypes.c_double, # mu  (<= 0 -> default 3.0e10 Pa)
    ctypes.c_int,    # na  (<= 0 -> default 5000)
    ctypes.c_double, # rmin (<= 0 -> default 2.0)
    ctypes.c_double, # rmax (<= 0 -> default 0.35)
    ctypes.c_int,    # n_surface (0 if no surface rupture)
    _int_p,          # surface_nodes[n_surface] 1-based
    ctypes.c_int,    # has_velocity (0 or 1)
    _dbl_p,          # velocity[ncells]
    ctypes.c_int,    # has_time (0 or 1)
    _dbl_p,          # time[nnodes] inout: initial conditions in, computed time out
    ctypes.c_int,    # has_pdf (0 or 1)
    _dbl_p,          # pdf_in[ncells]
    _dbl_p,          # slip[ncells] out
    _int_p,          # compute_time_flag out: 1 if time computed, 0 if slip only
]
_lib.k223d_compute_source.restype = None

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _dbl(arr):
    """Contiguous float64 array as a ctypes double pointer."""
    return arr.ctypes.data_as(_dbl_p)


def _int(arr):
    """Contiguous int32 array as a ctypes int pointer."""
    return arr.ctypes.data_as(_int_p)


def _dummy_dbl():
    """Single-element dummy double array for unused optional arguments."""
    return np.zeros(1, dtype=np.float64)


def _dummy_int():
    """Single-element dummy int32 array for unused optional arguments."""
    return np.zeros(1, dtype=np.int32)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_source(nodes, cells, mw,
                   mu=None, na=None, rmin=None, rmax=None,
                   velocity=None, time=None,
                   surface_nodes=None, pdf=None):
    """
    Compute a stochastic slip distribution on a triangular fault mesh.

    Parameters
    ----------
    nodes : array_like, shape (nnodes, 3)
        Node coordinates (x, y, z) in metres.
    cells : array_like, shape (ncells, 3), int
        Cell vertex indices, 1-based.
    mw : float
        Moment magnitude (required).
    mu : float, optional
        Shear modulus in Pa. Default: 3.0e10 Pa (30 GPa).
    na : int, optional
        Number of sub-events in the fractal cascade. Default: 5000.
    rmin : float, optional
        Minimum sub-event radius as a multiple of the mean cell radius.
        Default: 2.0. Should be >= 2.0 to resolve sub-events on the mesh.
    rmax : float, optional
        Maximum sub-event radius as a fraction of fault width.
        Default: 0.35. Must be < 0.5.
    velocity : array_like, shape (ncells,), optional
        Rupture velocity per cell in m/s. If provided together with `time`,
        rupture time is computed in addition to slip.
    time : array_like, shape (nnodes,), optional
        Initial rupture time per node in seconds.
        Convention: 0.0 at the nucleation node, 1.e32 at all other nodes.
        Updated in-place with computed rupture times when returned.
    surface_nodes : array_like of int, optional
        1-based indices of free-surface nodes. When provided, surface rupture
        is modelled with a reflected slip contribution.

    Returns
    -------
    slip : ndarray, shape (ncells,)
        Slip per cell in metres, normalised to the target seismic moment.
    time : ndarray, shape (nnodes,) or None
        Computed rupture time per node in seconds, or None if `velocity`
        and `time` were not both provided.

    Raises
    ------
    ValueError
        If array shapes are inconsistent.
    ImportError
        If libk223d.so cannot be found (raised at import time).

    Examples
    --------
    >>> slip, _ = k223d.compute_source(nodes, cells, mw=7.0)
    >>> slip, time = k223d.compute_source(nodes, cells, mw=7.0,
    ...                                    velocity=vel, time=t0)
    """
    # --- Validate and convert mesh arrays ---
    nodes = np.asarray(nodes, dtype=np.float64)
    cells = np.asarray(cells, dtype=np.int32)

    if nodes.ndim != 2 or nodes.shape[1] != 3:
        raise ValueError(f"nodes must have shape (nnodes, 3), got {nodes.shape}")
    if cells.ndim != 2 or cells.shape[1] != 3:
        raise ValueError(f"cells must have shape (ncells, 3), got {cells.shape}")

    if cells.min() == 0:
        # warnings.warn(
        #     "cells contains 0-based indices; converting to 1-based automatically. "
        #     "Pass cells+1 to silence this warning.",
        #     UserWarning, stacklevel=2
        # )
        cells = cells + 1

    nnodes = nodes.shape[0]
    ncells = cells.shape[0]

    # Ensure contiguous C-order layout before passing pointers
    nodes = np.ascontiguousarray(nodes)
    cells = np.ascontiguousarray(cells)

    # Extract coordinate columns as separate contiguous arrays
    px = np.ascontiguousarray(nodes[:, 0])
    py = np.ascontiguousarray(nodes[:, 1])
    pz = np.ascontiguousarray(nodes[:, 2])

    # --- Optional scalar parameters ---
    # Sentinel -1 tells Fortran to use the built-in default
    c_mu   = ctypes.c_double(mu   if mu   is not None else -1.0)
    c_na   = ctypes.c_int(   na   if na   is not None else -1)
    c_rmin = ctypes.c_double(rmin if rmin is not None else -1.0)
    c_rmax = ctypes.c_double(rmax if rmax is not None else -1.0)

    # --- Surface nodes ---
    if surface_nodes is not None:
        surf_arr = np.asarray(surface_nodes, dtype=np.int32)
        if len(surf_arr) == nnodes:
            # Flag array (one entry per node, 1=surface 0=interior): convert to index array
            surf_arr = (np.where(surf_arr > 0)[0] + 1).astype(np.int32)
        surf_arr = np.ascontiguousarray(surf_arr)
        n_surface = len(surf_arr)
    else:
        surf_arr = _dummy_int()
        n_surface = 0

    # --- Velocity ---
    if velocity is not None:
        vel_arr = np.ascontiguousarray(velocity, dtype=np.float64)
        if vel_arr.shape != (ncells,):
            raise ValueError(
                f"velocity must have shape ({ncells},), got {vel_arr.shape}"
            )
        has_vel = 1
    else:
        vel_arr = _dummy_dbl()
        has_vel = 0

    # --- Time (inout: initial conditions in, computed rupture time out) ---
    if time is not None:
        time_arr = np.ascontiguousarray(time, dtype=np.float64)
        if time_arr.shape != (nnodes,):
            raise ValueError(
                f"time must have shape ({nnodes},), got {time_arr.shape}"
            )
        has_time = 1
    else:
        time_arr = np.zeros(nnodes, dtype=np.float64)
        has_time = 0

    # --- PDF ---
    if pdf is not None:
        pdf_arr = np.ascontiguousarray(pdf, dtype=np.float64)
        if pdf_arr.shape != (ncells,):
            raise ValueError(
                f"pdf must have shape ({ncells},), got {pdf_arr.shape}"
            )
        has_pdf = 1
    else:
        pdf_arr = _dummy_dbl()
        has_pdf = 0

    # --- Output array ---
    slip_arr = np.zeros(ncells, dtype=np.float64)
    compute_time_flag = ctypes.c_int(0)

    # --- Call Fortran ---
    _lib.k223d_compute_source(
        ctypes.c_int(nnodes),
        ctypes.c_int(ncells),
        _dbl(px), _dbl(py), _dbl(pz),
        _int(cells),
        ctypes.c_double(float(mw)),
        c_mu, c_na, c_rmin, c_rmax,
        ctypes.c_int(n_surface),
        _int(surf_arr),
        ctypes.c_int(has_vel),
        _dbl(vel_arr),
        ctypes.c_int(has_time),
        _dbl(time_arr),
        ctypes.c_int(has_pdf),
        _dbl(pdf_arr),
        _dbl(slip_arr),
        ctypes.byref(compute_time_flag),
    )

    if compute_time_flag.value == 1:
        return slip_arr, time_arr
    else:
        return slip_arr, None
