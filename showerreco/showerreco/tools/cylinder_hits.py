"""
  showerreco.tools.cylinder_hit_finder
  -------------------------------------
  Python equivalent of Shower3DCylinderTrackHitFinder.

  Selects space points inside a cylinder of radius `max_perp_dist` and
  half-length `max_proj_dist` centred on the shower axis, starting from
  the shower start position.
"""

from __future__ import annotations

import numpy as np
from showerreco.geometry import ShowerElementHolder, SpacePoint

START_IN  = "ShowerStartPosition"
DIR_IN    = "ShowerDirection"
TRACK_SP_OUT   = "InitialTrackSpacePoints"
TRACK_IDX_OUT  = "InitialTrackIndices"   # indices into original spacepoint list

def run(
    spacepoints: list[SpacePoint],
    holder: ShowerElementHolder,
    *,
    max_proj_dist:  float = 10.0,    # cm along axis  (fMaxProjectionDist)
    max_perp_dist:  float = 1.0,    # cm transverse  (fMaxPerpendicularDist)
    forward_only:   bool  = True,   # only downstream hits (fForwardHitsOnly)
    verbose: bool = False,
) -> int:

    if not holder.check_element(START_IN):
        if verbose: print("[CylinderHitFinder] Start position not set.")
        return 1
    if not holder.check_element(DIR_IN):
        if verbose: print("[CylinderHitFinder] Direction not set.")
        return 1
    if len(spacepoints) == 0:
        if verbose: print("[CylinderHitFinder] No spacepoints.")
        return 1

    start = holder.get_element(START_IN)
    dirn  = holder.get_element(DIR_IN)
    u     = dirn / np.linalg.norm(dirn)

    positions = np.array([sp.position for sp in spacepoints])
    delta     = positions - start

    proj = delta @ u                                          # longitudinal
    perp = np.linalg.norm(delta - np.outer(proj, u), axis=1) # transverse

    mask = (np.abs(proj) < max_proj_dist) & (np.abs(perp) < max_perp_dist)
    if forward_only:
        mask &= (proj >= 0)

    indices = np.where(mask)[0]
    selected = [spacepoints[i] for i in indices]

    holder.set_element(TRACK_SP_OUT,  selected)
    holder.set_element(TRACK_IDX_OUT, indices)

    if verbose:
        print(f"[CylinderHitFinder] Selected {len(selected)} / {len(spacepoints)} spacepoints "
              f"(proj<{max_proj_dist} cm, perp<{max_perp_dist} cm)")
    return 0