"""
    showerreco.tools.sliding_fit
    -----------------------------
    Python equivalent of ShowerPandoraSlidingFitTrackFinder.

    ICARUS configuration (showerpandoraslidingfittrackfinder):
        SlidingFitHalfWindow : 12 wire pitches × 0.3 cm = 3.6 cm physical half-window
        MinTrajectoryPoints  : 2
        dEdxTrackLength      : 3 cm  (only points within this range matter downstream)

    Since the half-window (3.6 cm) exceeds the dEdx stem length (3 cm), every
    window contains essentially all cylinder points — the fit is effectively a
    single global straight line anchored at the start, not a truly local sliding
    fit.  We implement this faithfully: the window is defined in physical cm so
    the behaviour is transparent and the trajectory points are the sorted cylinder
    spacepoints each assigned the same (global) fitted direction.

    Output stored in the holder under "SlidingFitTrack":
        dict with keys:
            'positions'   np.ndarray (N, 3)  — one per cylinder spacepoint, sorted by proj
            'directions'  np.ndarray (N, 3)  — local unit direction at each point
            'proj'        np.ndarray (N,)    — signed distance from start along PCA axis
            'fit_line'    dict with 'point' and 'direction' — the global fit for plotting
"""

from __future__ import annotations
import numpy as np
from showerreco.geometry import ShowerElementHolder


START_IN      = "ShowerStartPosition"
DIR_IN        = "ShowerDirection"
CYLINDER_IN   = "InitialTrackSpacePoints"
TRACK_OUT     = "SlidingFitTrack"
TRACK_LEN_OUT = "SlidingFitTrackLength"


def run(
    holder: ShowerElementHolder,
    *,
    half_window_cm:  float = 3.6,   # SlidingFitHalfWindow=12 × wire_pitch=0.3 cm
    dedx_length_cm:  float = 3.0,   # dEdxTrackLength: only keep points within this
    min_traj_points: int   = 2,     # MinTrajectoryPoints
    wire_pitch:      float = 0.3,   # ICARUS W-plane wire pitch in cm
    verbose: bool = False,
) -> int:
    """
        Parameters
        ----------
        half_window_cm  : Physical half-width of the sliding window in cm.
                        Set from SlidingFitHalfWindow × wire_pitch.
        dedx_length_cm  : Maximum longitudinal distance from start to include a
                        trajectory point (mirrors dEdxTrackLength).
        min_traj_points : Minimum trajectory points to declare success.
        wire_pitch      : Wire pitch used to convert SlidingFitHalfWindow to cm.
    """

    for key in (START_IN, DIR_IN, CYLINDER_IN):
        if not holder.check_element(key):
            if verbose:
                print(f"[SlidingFit] '{key}' not set — returning.")
            return 1

    start = holder.get_element(START_IN)
    dirn  = holder.get_element(DIR_IN)
    u     = dirn / np.linalg.norm(dirn)

    cylinder_sps = holder.get_element(CYLINDER_IN)
    if len(cylinder_sps) < min_traj_points:
        if verbose:
            print(f"[SlidingFit] Only {len(cylinder_sps)} spacepoints — insufficient.")
        return 1

    positions = np.array([sp.position for sp in cylinder_sps])  # (N, 3)

    # --- sort by longitudinal projection from start -----------------------
    proj_vals = (positions - start) @ u
    order     = np.argsort(proj_vals)
    positions = positions[order]
    proj_vals = proj_vals[order]

    # --- global linear fit over all cylinder points -----------------------
    # (since half_window >= stem length this is what every window reduces to)
    centre_global = positions.mean(axis=0)
    delta_global  = positions - centre_global
    cov_global    = (delta_global.T @ delta_global) / len(positions)
    evals, evecs  = np.linalg.eigh(cov_global)
    global_dir    = evecs[:, -1]
    if np.dot(global_dir, u) < 0:
        global_dir = -global_dir

    fit_line = {"point": centre_global, "direction": global_dir}

    # --- sliding window: cm-based, one trajectory point per spacepoint ----
    N = len(positions)
    traj_positions  = []
    traj_directions = []
    traj_proj       = []

    for i in range(N):
        s_i = proj_vals[i]

        # skip points beyond the dEdx stem length
        if s_i > dedx_length_cm:
            continue

        # select window by physical distance
        window_mask = np.abs(proj_vals - s_i) <= half_window_cm
        window      = positions[window_mask]

        if len(window) < 2:
            # fall back to global direction
            local_dir = global_dir.copy()
            centre    = positions[i]
        else:
            centre   = window.mean(axis=0)
            delta_w  = window - centre
            cov_w    = (delta_w.T @ delta_w) / len(window)
            ev, evec = np.linalg.eigh(cov_w)
            local_dir = evec[:, -1]
            if np.dot(local_dir, u) < 0:
                local_dir = -local_dir

        traj_positions.append(positions[i])   # position = the actual spacepoint
        traj_directions.append(local_dir)
        traj_proj.append(s_i)

    if len(traj_positions) < min_traj_points:
        if verbose:
            print(f"[SlidingFit] Only {len(traj_positions)} trajectory points — insufficient.")
        return 1

    traj_positions  = np.array(traj_positions)
    traj_directions = np.array(traj_directions)
    traj_proj       = np.array(traj_proj)

    track = {
        "positions":  traj_positions,
        "directions": traj_directions,
        "proj":       traj_proj,
        "fit_line":   fit_line,        # global fit, useful for plotting
    }
    holder.set_element(TRACK_OUT, track)

    track_length = float(traj_proj[-1] - traj_proj[0]) if len(traj_proj) > 1 else 0.
    holder.set_element(TRACK_LEN_OUT, track_length)

    if verbose:
        angle_spread = np.degrees(np.arccos(np.clip(
            [np.dot(d, global_dir) for d in traj_directions], -1, 1
        )))
        print(f"[SlidingFit] {len(traj_positions)} trajectory points within "
              f"{dedx_length_cm} cm stem")
        print(f"[SlidingFit] Track length = {track_length:.2f} cm")
        print(f"[SlidingFit] Half-window  = {half_window_cm:.1f} cm  "
              f"({half_window_cm/wire_pitch:.0f} wire pitches)")
        print(f"[SlidingFit] Local direction spread: "
              f"max {angle_spread.max():.2f}°, mean {angle_spread.mean():.2f}°")

    return 0