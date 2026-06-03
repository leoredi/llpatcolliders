"""
event_display.py

Three-panel event displays (3D + top-down x-z + side y-z) for GRENDEL
signal MC events. Driven from decayProbPerEvent_2body.py: pass the MC
output dict together with a dictionary of named per-sample selections,
and the first N events satisfying each named selection are drawn,
picking one representative MC sample per particle (highest-weight
sample inside the selection mask).

API
---
    make_event_displays(mc, selections, n_per=3,
                        out_prefix='event_display', out_dir='.')

Each event shows:
  - the cavern centreline + a handful of cross-section ribs (3D)
  - the IP origin
  - the LLP trajectory from the IP to the decay vertex
  - the decay vertex marker
  - both daughter rays from the vertex to their layer-1 wall hit
  - the predicted layer-2 hit (extension by L = DETECTOR_THICKNESS
    radially) shown as a dashed segment when defined
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers '3d')

from grendel_geometry import (
    DETECTOR_THICKNESS, path_3d_fiducial,
    classify_points_with_basis,
    _profile_pts as profile_pts,
)


_DAUGHTER_COLORS = ('tab:blue', 'tab:green')
_PARTICLE_COLORS = ('tab:orange', 'tab:purple', 'tab:red', 'tab:brown')


def _local_basis(i, p):
    """Local tunnel basis at centreline index i (tangent, right, up)."""
    if i == 0:
        tan = p[1] - p[0]
    elif i == len(p) - 1:
        tan = p[-1] - p[-2]
    else:
        tan = p[i + 1] - p[i - 1]
    tan = tan / np.linalg.norm(tan)
    wu = np.array([0., 1., 0.]) if abs(tan[1]) < 0.9 \
        else np.array([0., 0., 1.])
    right = np.cross(tan, wu); right /= np.linalg.norm(right)
    up = np.cross(right, tan); up /= np.linalg.norm(up)
    return tan, right, up


def _project_3d(ax3d, pts3, **kwargs):
    pts3 = np.atleast_2d(pts3)
    ax3d.plot(pts3[:, 0], pts3[:, 2], pts3[:, 1], **kwargs)


def _scatter_3d(ax3d, p, **kwargs):
    ax3d.scatter([p[0]], [p[2]], [p[1]], **kwargs)


def _draw_cavern(ax3d, ax_xz, ax_yz, n_ribs=8):
    """Centreline in all three views + cross-section ribs in 3D + IP."""
    p = path_3d_fiducial

    # Centrelines
    _project_3d(ax3d, p, color='gray', alpha=0.5, lw=1.2)
    ax_xz.plot(p[:, 2], p[:, 0], color='gray', alpha=0.5, lw=1.2)
    ax_yz.plot(p[:, 2], p[:, 1], color='gray', alpha=0.5, lw=1.2)

    # IP origin
    ax3d.scatter([0], [0], [0], color='crimson', marker='*', s=70,
                 label='IP', zorder=10)
    ax_xz.scatter([0], [0], color='crimson', marker='*', s=70,
                  zorder=10, label='IP')
    ax_yz.scatter([0], [0], color='crimson', marker='*', s=70,
                  zorder=10, label='IP')

    # 3D cross-section ribs
    idxs = np.linspace(0, len(p) - 1, n_ribs).astype(int)
    for i in idxs:
        _, right, up = _local_basis(i, p)
        ring = p[i] + np.outer(profile_pts[:, 0], right) \
            + np.outer(profile_pts[:, 1], up)
        ring = np.vstack([ring, ring[0:1]])
        _project_3d(ax3d, ring, color='gray', alpha=0.3, lw=0.5)

    # 2D envelopes — for each centreline point, project the profile
    # outline. Build a polygon trace per view by collecting min/max in
    # the projected coordinate at each centreline point.
    x_low, x_high = [], []
    y_low, y_high = [], []
    z_centre = []
    for i in range(len(p)):
        _, right, up = _local_basis(i, p)
        ring = p[i] + np.outer(profile_pts[:, 0], right) \
            + np.outer(profile_pts[:, 1], up)
        x_low.append(ring[:, 0].min())
        x_high.append(ring[:, 0].max())
        y_low.append(ring[:, 1].min())
        y_high.append(ring[:, 1].max())
        z_centre.append(p[i, 2])
    z_centre = np.array(z_centre)
    ax_xz.fill_between(z_centre, x_low, x_high, color='gray', alpha=0.10)
    ax_yz.fill_between(z_centre, y_low, y_high, color='gray', alpha=0.10)


def _draw_particle(mc, idx, ax3d, ax_xz, ax_yz, color):
    """Draw one particle (one representative sample) into the 3 axes."""
    decay_pos = mc['decay_pos'][idx]
    dir1 = mc['dir1'][idx]
    dir2 = mc['dir2'][idx]
    p1_1 = mc['exit_pt_1'][idx]
    p1_2 = mc['exit_pt_2'][idx]

    # LLP trajectory: IP -> decay vertex
    llp = np.vstack([np.zeros(3), decay_pos])
    _project_3d(ax3d, llp, color=color, lw=1.0, alpha=0.7)
    ax_xz.plot(llp[:, 2], llp[:, 0], color=color, lw=1.0, alpha=0.7)
    ax_yz.plot(llp[:, 2], llp[:, 1], color=color, lw=1.0, alpha=0.7)

    # Decay vertex
    _scatter_3d(ax3d, decay_pos, color=color, marker='o', s=35, zorder=8)
    ax_xz.scatter([decay_pos[2]], [decay_pos[0]],
                  color=color, marker='o', s=35, zorder=8)
    ax_yz.scatter([decay_pos[2]], [decay_pos[1]],
                  color=color, marker='o', s=35, zorder=8)

    # Daughter rays + layer-1/layer-2 hits
    for (p_hit, d_vec, dcolor) in ((p1_1, dir1, _DAUGHTER_COLORS[0]),
                                    (p1_2, dir2, _DAUGHTER_COLORS[1])):
        if not np.all(np.isfinite(p_hit)):
            continue
        # vertex -> layer-1 hit
        seg1 = np.vstack([decay_pos, p_hit])
        _project_3d(ax3d, seg1, color=dcolor, lw=1.1)
        ax_xz.plot(seg1[:, 2], seg1[:, 0], color=dcolor, lw=1.1)
        ax_yz.plot(seg1[:, 2], seg1[:, 1], color=dcolor, lw=1.1)
        _scatter_3d(ax3d, p_hit, color=dcolor, marker='s', s=30, zorder=8)
        ax_xz.scatter([p_hit[2]], [p_hit[0]],
                      color=dcolor, marker='s', s=30, zorder=8)
        ax_yz.scatter([p_hit[2]], [p_hit[1]],
                      color=dcolor, marker='s', s=30, zorder=8)

        # Predicted layer-2 hit (parallel-plane extension)
        theta, _, _, right, up = classify_points_with_basis(p_hit[None, :])
        n_hat = np.cos(theta[0]) * right[0] + np.sin(theta[0]) * up[0]
        d_dot_n = float(np.dot(d_vec, n_hat))
        if abs(d_dot_n) > 5e-2:
            p2 = p_hit + (DETECTOR_THICKNESS / d_dot_n) * d_vec
            seg2 = np.vstack([p_hit, p2])
            _project_3d(ax3d, seg2, color=dcolor, lw=1.1, linestyle='--')
            ax_xz.plot(seg2[:, 2], seg2[:, 0],
                       color=dcolor, lw=1.1, linestyle='--')
            ax_yz.plot(seg2[:, 2], seg2[:, 1],
                       color=dcolor, lw=1.1, linestyle='--')
            _scatter_3d(ax3d, p2, color=dcolor, marker='D', s=25, zorder=8)
            ax_xz.scatter([p2[2]], [p2[0]],
                          color=dcolor, marker='D', s=25, zorder=8)
            ax_yz.scatter([p2[2]], [p2[1]],
                          color=dcolor, marker='D', s=25, zorder=8)


def make_event_displays(mc, selections, n_per=3,
                        out_prefix='event_display', out_dir='.'):
    """
    Save 3-panel event displays for the first ``n_per`` events satisfying
    each named selection.

    Parameters
    ----------
    mc : dict
        sample_separations() output. Must include the per-(sample) arrays
        ``decay_pos``, ``dir1``, ``dir2``, ``exit_pt_1``, ``exit_pt_2``,
        ``event``, ``pid``, ``weights``, plus the cut variables.
    selections : dict[str, np.ndarray[bool]]
        Per-sample boolean masks. The dict keys become filename suffixes
        and are used to label each figure. Edit this dict in the caller
        to control which event populations get displayed.
    n_per : int
        First N events per selection.
    out_prefix : str
        Output filename prefix.
    out_dir : str
        Output directory.
    """
    if len(mc.get('decay_pos', [])) == 0:
        print("[event_display] empty MC; nothing to draw")
        return

    pid = mc['pid']
    event = mc['event']
    weights = mc['weights']

    os.makedirs(out_dir, exist_ok=True)

    for sel_name, sel_mask in selections.items():
        sel_mask = np.asarray(sel_mask, dtype=bool)
        if not sel_mask.any():
            print(f"[event_display] '{sel_name}' has no passing samples; "
                  f"skip")
            continue

        ev_pass = np.unique(event[sel_mask])
        ev_pass.sort()
        ev_show = ev_pass[:n_per]
        print(f"[event_display] '{sel_name}': drawing {len(ev_show)} "
              f"events (of {len(ev_pass)} candidates)")

        for k, ev in enumerate(ev_show):
            ev_mask = (event == ev)
            particles = np.unique(pid[ev_mask])

            # Pick one representative sample per particle inside the
            # selection mask
            rep_idxs = []
            for p_id in particles:
                p_mask = (pid == p_id) & ev_mask & sel_mask
                if not p_mask.any():
                    continue
                rep_idxs.append(int(np.argmax(np.where(p_mask, weights, 0))))
            if not rep_idxs:
                continue

            fig = plt.figure(figsize=(18, 6.5))
            ax3d = fig.add_subplot(1, 3, 1, projection='3d')
            ax_xz = fig.add_subplot(1, 3, 2)
            ax_yz = fig.add_subplot(1, 3, 3)

            _draw_cavern(ax3d, ax_xz, ax_yz)
            for j, idx in enumerate(rep_idxs):
                _draw_particle(mc, idx, ax3d, ax_xz, ax_yz,
                               color=_PARTICLE_COLORS[j % len(_PARTICLE_COLORS)])

            ax3d.set_xlabel('x (m)')
            ax3d.set_ylabel('z (m)')
            ax3d.set_zlabel('y (m)')
            ax3d.set_title('3D')

            ax_xz.set_xlabel('z (m)')
            ax_xz.set_ylabel('x (m)')
            ax_xz.set_aspect('equal')
            ax_xz.grid(True, alpha=0.3)
            ax_xz.set_title('Top-down (x vs z)')

            ax_yz.set_xlabel('z (m)')
            ax_yz.set_ylabel('y (m)')
            ax_yz.set_aspect('equal')
            ax_yz.grid(True, alpha=0.3)
            ax_yz.set_title('Side (y vs z)')

            header = [
                f"Event {int(ev)}    selection: {sel_name}    "
                f"({k + 1}/{len(ev_show)})"
            ]
            for j, idx in enumerate(rep_idxs):
                header.append(
                    f"  p{j}: d={mc['d'][idx]:.2f} m, "
                    f"sep_in={mc['sep'][idx]*100:.1f} cm, "
                    f"sep_out={mc['sep_outer'][idx]*100:.1f} cm, "
                    f"pointing={mc['pointing'][idx]*1000:.0f} mrad, "
                    f"dca={mc['dca'][idx]*100:.2f} cm"
                )
            fig.suptitle("\n".join(header), fontsize=9, x=0.02, ha='left')
            fig.tight_layout(rect=(0, 0, 1, 0.93))

            fname = os.path.join(
                out_dir, f"{out_prefix}_{sel_name}_{k:02d}.png")
            fig.savefig(fname, dpi=120)
            plt.close(fig)
            print(f"  saved {fname}")
