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
    DETECTOR_THICKNESS, path_3d_fiducial, cumulative_length,
    classify_points_with_basis, tunnel_profile_points,
    _profile_pts as profile_pts,
)

# Inner-tracker-surface (= fiducial mesh boundary) and outer
# (= original tunnel wall) profiles, in the same local frame.
_profile_layer1 = tunnel_profile_points(inset=DETECTOR_THICKNESS)
_profile_layer2 = profile_pts


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


def _draw_cavern(ax3d, ax_xz, ax_yz, bounds, n_ribs=6):
    """Centreline + cross-section ribs + outline within ``bounds`` only."""
    p = path_3d_fiducial
    xlim, ylim, zlim = bounds

    in_box = ((p[:, 0] >= xlim[0]) & (p[:, 0] <= xlim[1]) &
              (p[:, 2] >= zlim[0]) & (p[:, 2] <= zlim[1]))
    idx = np.where(in_box)[0]
    if len(idx) == 0:
        return
    lo = max(0, int(idx[0]) - 1)
    hi = min(len(p), int(idx[-1]) + 2)
    p_local = p[lo:hi]

    # Centrelines
    _project_3d(ax3d, p_local, color='gray', alpha=0.5, lw=1.2)
    ax_xz.plot(p_local[:, 2], p_local[:, 0], color='gray', alpha=0.5, lw=1.2)
    ax_yz.plot(p_local[:, 2], p_local[:, 1], color='gray', alpha=0.5, lw=1.2)

    # 3D ribs at evenly-spaced indices in the local subset
    rib_ix = np.linspace(0, len(p_local) - 1,
                          min(n_ribs, len(p_local))).astype(int)
    for j in rib_ix:
        gi = lo + int(j)
        _, right, up = _local_basis(gi, p)
        ring = p[gi] + np.outer(profile_pts[:, 0], right) \
            + np.outer(profile_pts[:, 1], up)
        ring = np.vstack([ring, ring[0:1]])
        _project_3d(ax3d, ring, color='gray', alpha=0.3, lw=0.5)

    # 2D envelopes over the local subset only
    x_low, x_high, y_low, y_high, z_centre = [], [], [], [], []
    for i in range(lo, hi):
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


def _zoom_bounds(mc, idxs, padding=3.0):
    """3D bounding box (xlim, ylim, zlim) around an event's key points."""
    pts = []
    for idx in idxs:
        pts.append(mc['decay_pos'][idx])
        for k in (1, 2):
            p = mc[f'exit_pt_{k}'][idx]
            if np.all(np.isfinite(p)):
                pts.append(p)
    pts = np.array(pts)
    return (
        (float(pts[:, 0].min() - padding), float(pts[:, 0].max() + padding)),
        (float(pts[:, 1].min() - padding), float(pts[:, 1].max() + padding)),
        (float(pts[:, 2].min() - padding), float(pts[:, 2].max() + padding)),
    )


def _centerline_at(s):
    """Centreline point + basis at arc-length s along the fiducial path."""
    seg_idx = int(np.searchsorted(cumulative_length, s) - 1)
    seg_idx = max(0, min(seg_idx, len(path_3d_fiducial) - 2))
    seg = path_3d_fiducial[seg_idx + 1] - path_3d_fiducial[seg_idx]
    seg_len = np.linalg.norm(seg)
    seg_hat = seg / seg_len
    t = max(0.0, min(seg_len, s - cumulative_length[seg_idx]))
    cl_point = path_3d_fiducial[seg_idx] + t * seg_hat
    wu = np.array([0., 1., 0.]) if abs(seg_hat[1]) < 0.9 \
        else np.array([0., 0., 1.])
    right = np.cross(seg_hat, wu); right /= np.linalg.norm(right)
    up = np.cross(right, seg_hat); up /= np.linalg.norm(up)
    return cl_point, seg_hat, right, up


def _draw_local_section(mc, idxs, ax_local):
    """Bottom panel: cavern cross-section + layer 1/2 at the wall hits.

    The cross-section is drawn at the average DECAY VERTEX's s along the
    centreline (not the wall hit's s). This way the decay vertex sits in
    its own perpendicular plane and lands inside the fiducial outline as
    it physically should. The wall hits are typically at slightly larger
    s (daughters travel forward), and are projected onto the decay's
    cross-section for context — for long LLP paths through the fiducial
    they may not coincide exactly with the drawn profile.
    """
    hits3d = []
    for idx in idxs:
        for k in (1, 2):
            p = mc[f'exit_pt_{k}'][idx]
            if np.all(np.isfinite(p)):
                hits3d.append(p)
    if not hits3d:
        ax_local.set_title('Local cross-section (no wall hits)')
        return
    hits3d = np.array(hits3d)

    # Reference s = average decay vertex along the centreline, so decays
    # are drawn in their own cross-section.
    decay_arr = np.array([mc['decay_pos'][idx] for idx in idxs])
    decay_avg = decay_arr.mean(axis=0)
    _, s_dec, _, _, _ = classify_points_with_basis(decay_avg[None, :])
    s_ref = float(s_dec[0])

    # For context, also report the average wall-hit s.
    hit_avg = hits3d.mean(axis=0)
    _, s_hit_arr, _, _, _ = classify_points_with_basis(hit_avg[None, :])
    s_hit_avg = float(s_hit_arr[0])

    cl_point, _, right, up = _centerline_at(s_ref)

    def to_local(p3d):
        rel = p3d - cl_point
        return np.array([np.dot(rel, right), np.dot(rel, up)])

    # Draw the two profile contours (layer 1 inner = fiducial boundary,
    # layer 2 outer = tunnel wall at +L radial)
    ring1 = np.vstack([_profile_layer1, _profile_layer1[0:1]])
    ring2 = np.vstack([_profile_layer2, _profile_layer2[0:1]])
    ax_local.plot(ring1[:, 0], ring1[:, 1], color='black', lw=1.5,
                  label='Layer 1 (fiducial wall)')
    ax_local.plot(ring2[:, 0], ring2[:, 1], color='gray', lw=1.5,
                  linestyle='--', label=f'Layer 2 (+{DETECTOR_THICKNESS*100:.0f} cm)')

    # Per-particle: decay vertex (star), LLP back-trace, both daughter
    # rays from the vertex out through layer-1 then on to layer-2.
    extra_xys = []
    for j, idx in enumerate(idxs):
        decay_pos = mc['decay_pos'][idx]
        decay_xy = to_local(decay_pos)
        pcolor = _PARTICLE_COLORS[j % len(_PARTICLE_COLORS)]

        # LLP back-trace: line from a point 3 m behind the vertex along
        # -LLP direction up to the vertex.
        llp_norm = np.linalg.norm(decay_pos)
        if llp_norm > 1e-6:
            llp_dir3d = decay_pos / llp_norm
            llp_dir2d = np.array([np.dot(llp_dir3d, right),
                                   np.dot(llp_dir3d, up)])
            llp_back = decay_xy - 3.0 * llp_dir2d
            ax_local.plot([llp_back[0], decay_xy[0]],
                          [llp_back[1], decay_xy[1]],
                          color=pcolor, lw=1.0, alpha=0.7,
                          linestyle=':')
            extra_xys.append(llp_back)

        # Decay vertex marker
        ax_local.scatter(decay_xy[0], decay_xy[1],
                         color=pcolor, marker='*', s=90, zorder=10,
                         edgecolors='black', linewidths=0.5)
        extra_xys.append(decay_xy)

        for k, dcolor in ((1, _DAUGHTER_COLORS[0]),
                          (2, _DAUGHTER_COLORS[1])):
            p_hit = mc[f'exit_pt_{k}'][idx]
            d_vec = mc[f'dir{k}'][idx]
            if not np.all(np.isfinite(p_hit)):
                continue
            p_hit_xy = to_local(p_hit)
            d_xy = np.array([np.dot(d_vec, right), np.dot(d_vec, up)])

            # Daughter trajectory from vertex to layer-1 hit
            ax_local.plot([decay_xy[0], p_hit_xy[0]],
                          [decay_xy[1], p_hit_xy[1]],
                          color=dcolor, lw=1.2)
            ax_local.scatter(p_hit_xy[0], p_hit_xy[1],
                             color=dcolor, marker='s', s=45, zorder=8)

            # Predicted layer-2 hit
            theta_p = np.arctan2(p_hit_xy[1], p_hit_xy[0])
            n_hat_2d = np.array([np.cos(theta_p), np.sin(theta_p)])
            d_dot_n = float(np.dot(d_xy, n_hat_2d))
            if abs(d_dot_n) > 5e-2:
                p2_xy = p_hit_xy + (DETECTOR_THICKNESS / d_dot_n) * d_xy
                ax_local.plot([p_hit_xy[0], p2_xy[0]],
                              [p_hit_xy[1], p2_xy[1]],
                              color=dcolor, lw=1.2, linestyle='--')
                ax_local.scatter(p2_xy[0], p2_xy[1],
                                 color=dcolor, marker='D', s=30, zorder=8)

    # Zoom to include hits, decay vertex, and LLP back-trace start
    all_xys = np.vstack([np.array([to_local(p) for p in hits3d]),
                          np.array(extra_xys)])
    cx, cy = all_xys.mean(axis=0)
    spread = max(float(np.ptp(all_xys[:, 0])),
                 float(np.ptp(all_xys[:, 1])))
    span = max(0.6, 0.6 + spread / 2)
    ax_local.set_xlim(cx - span, cx + span)
    ax_local.set_ylim(cy - span, cy + span)
    ax_local.set_aspect('equal')
    ax_local.set_xlabel('x_local (m)')
    ax_local.set_ylabel('y_local (m)')
    ds = s_hit_avg - s_ref
    ax_local.set_title(
        f'Local cross-section at decay s = {s_ref:.1f} m '
        f'(hit s = {s_hit_avg:.1f} m, Δs = {ds:+.1f} m; zoom ±{span:.1f} m)')
    ax_local.legend(fontsize=8, loc='upper right')
    ax_local.grid(True, alpha=0.3)


def _draw_particle(mc, idx, ax3d, ax_xz, ax_yz, color):
    """Draw one particle (one representative sample) into the 3 axes."""
    decay_pos = mc['decay_pos'][idx]
    dir1 = mc['dir1'][idx]
    dir2 = mc['dir2'][idx]
    p1_1 = mc['exit_pt_1'][idx]
    p1_2 = mc['exit_pt_2'][idx]

    # LLP trajectory (last 3 m before the decay vertex — the IP itself
    # is outside the cavern-zoomed view).
    llp_norm = np.linalg.norm(decay_pos)
    if llp_norm > 1e-6:
        llp_dir = decay_pos / llp_norm
        llp_back = decay_pos - 3.0 * llp_dir
        llp = np.vstack([llp_back, decay_pos])
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
                        out_prefix='event_display', out_dir='.',
                        rng_seed=42):
    """
    Save 3-panel event displays for ``n_per`` events randomly drawn
    (without replacement) from each named selection.

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
        Number of events to sample per selection.
    out_prefix : str
        Output filename prefix.
    out_dir : str
        Output directory.
    rng_seed : int
        Seed for the random event picker — change to get a different draw.
    """
    if len(mc.get('decay_pos', [])) == 0:
        print("[event_display] empty MC; nothing to draw")
        return

    pid = mc['pid']
    event = mc['event']
    weights = mc['weights']
    rng = np.random.default_rng(rng_seed)

    os.makedirs(out_dir, exist_ok=True)

    for sel_name, sel_mask in selections.items():
        sel_mask = np.asarray(sel_mask, dtype=bool)
        if not sel_mask.any():
            print(f"[event_display] '{sel_name}' has no passing samples; "
                  f"skip")
            continue

        ev_pass = np.unique(event[sel_mask])
        n_draw = min(n_per, len(ev_pass))
        ev_show = rng.choice(ev_pass, size=n_draw, replace=False)
        ev_show.sort()
        print(f"[event_display] '{sel_name}': drawing {len(ev_show)} "
              f"random events (of {len(ev_pass)} candidates)")

        for k, ev in enumerate(ev_show):
            ev_mask = (event == ev)
            particles = np.unique(pid[ev_mask])

            # Pick one representative sample per particle inside the
            # selection mask. Sample weight-proportionally — picking the
            # argmax instead biases the representative to the smallest-d
            # passing sample (since w ∝ exp(-d/λ) is monotone in d, the
            # max sits at the lowest d, putting the displayed vertex
            # right at LLP entry for every event).
            rep_idxs = []
            for p_id in particles:
                p_mask = (pid == p_id) & ev_mask & sel_mask
                if not p_mask.any():
                    continue
                idxs_p = np.where(p_mask)[0]
                w_p = weights[idxs_p]
                if w_p.sum() <= 0:
                    rep_idxs.append(int(rng.choice(idxs_p)))
                else:
                    rep_idxs.append(int(rng.choice(idxs_p, p=w_p / w_p.sum())))
            if not rep_idxs:
                continue

            fig = plt.figure(figsize=(18, 12))
            gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0])
            ax3d = fig.add_subplot(gs[0, 0], projection='3d')
            ax_xz = fig.add_subplot(gs[0, 1])
            ax_yz = fig.add_subplot(gs[0, 2])
            ax_local = fig.add_subplot(gs[1, :])

            bounds = _zoom_bounds(mc, rep_idxs, padding=3.0)
            xlim, ylim, zlim = bounds

            _draw_cavern(ax3d, ax_xz, ax_yz, bounds)
            for j, idx in enumerate(rep_idxs):
                _draw_particle(mc, idx, ax3d, ax_xz, ax_yz,
                               color=_PARTICLE_COLORS[j % len(_PARTICLE_COLORS)])
            _draw_local_section(mc, rep_idxs, ax_local)

            ax3d.set_xlim(*xlim)
            ax3d.set_ylim(*zlim)
            ax3d.set_zlim(*ylim)
            ax3d.set_xlabel('x (m)')
            ax3d.set_ylabel('z (m)')
            ax3d.set_zlabel('y (m)')
            ax3d.set_title('3D (zoom)')

            ax_xz.set_aspect('auto')
            ax_xz.set_xlim(*zlim)
            ax_xz.set_ylim(*xlim)
            ax_xz.set_xlabel('z (m)')
            ax_xz.set_ylabel('x (m)')
            ax_xz.grid(True, alpha=0.3)
            ax_xz.set_title('Top-down (x vs z)')

            ax_yz.set_aspect('auto')
            ax_yz.set_xlim(*zlim)
            ax_yz.set_ylim(*ylim)
            ax_yz.set_xlabel('z (m)')
            ax_yz.set_ylabel('y (m)')
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
