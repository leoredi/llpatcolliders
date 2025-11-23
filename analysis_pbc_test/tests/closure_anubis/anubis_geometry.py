"""
ANUBIS vertical shaft geometry following Hirsch & Wang (2001.04750v2).

Reference: Fig. 1 of arXiv:2001.04750v2 [hep-ph]

- IP (pp collision point) at (0, 0, 0).
- Vertical shaft is a cylinder with:
    * radius  R = 9.0 m (diameter lh = 18 m)
    * height  H = 56.0 m (length lv = 56 m)
    * horizontal displacement from IP: dh = 5.0 m
    * vertical displacement from IP: dv = 24.0 m
    * centre  C = (5.0, 24.0 + 28.0, 0.0) = (5.0, 52.0, 0.0) m
  => bottom at y = 24.0 m, top at y = 80.0 m

Axis is parallel to +y (vertical).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import trimesh

# Ensure analysis_pbc_test is on sys.path so we can import geometry.per_parent_efficiency
THIS_FILE = Path(__file__).resolve()
# .../analysis_pbc_test/tests/closure_anubis/anubis_geometry.py
ANALYSIS_DIR = THIS_FILE.parents[2]  # .../analysis_pbc_test
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from geometry.per_parent_efficiency import create_tube_mesh


def build_anubis_shaft_mesh(
    radius_m: float = 9.0,        # lh/2 = 18/2 = 9 m
    height_m: float = 56.0,       # lv = 56 m
    center_x_m: float = 5.0,      # dh = 5 m (horizontal displacement)
    center_y_m: float = 52.0,     # dv + lv/2 = 24 + 28 = 52 m
    center_z_m: float = 0.0,
    n_segments: int = 32,
) -> trimesh.Trimesh:
    """
    Build a simplified ANUBIS-like vertical shaft as a trimesh.Trimesh.

    Parameters
    ----------
    radius_m : float
        Cylinder radius in metres.
    height_m : float
        Cylinder height in metres.
    center_x_m, center_y_m, center_z_m : float
        Coordinates of the cylinder centre in metres.
    n_segments : int
        Number of azimuthal segments for the tube mesh.

    Returns
    -------
    mesh : trimesh.Trimesh
        Volume mesh of the cylindrical shaft.
    """
    half_h = 0.5 * height_m

    # Bottom and top points along the +y axis
    bottom = np.array([center_x_m, center_y_m - half_h, center_z_m], dtype=float)
    top = np.array([center_x_m, center_y_m + half_h, center_z_m], dtype=float)

    path_points = np.stack([bottom, top], axis=0)

    # create_tube_mesh returns (vertices, faces)
    vertices, faces = create_tube_mesh(
        path_points=path_points,
        radius=radius_m,
        n_segments=n_segments,
    )

    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
    return mesh
