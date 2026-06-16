"""
    showerreco.tools.start_position
    -------------------------------
    Python equivalent of ShowerPFPVertexStartPosition.

    Original C++ logic (ShowerPFPVertexStartPosition.cc):
    1. Retrieve the recob::Vertex associated with the PFParticle.
    2. If exactly one vertex exists, store it as the shower start position.
    3. Fall back to ordering space-points along the PCA direction if no vertex
        is available (requires ShowerDirection to already be set).

    Here the "PFP vertex" is simply the neutrino interaction vertex that Pandora
    attaches to each shower PFParticle.  In our Python representation we receive
    it directly as a Point, which is the common use-case for event-display work.
"""

from __future__ import annotations

import numpy

from showerreco.geometry import (
    Point, Vector, SpacePoint,
    ShowerElementHolder,
    UNSET_POINT, UNSET_VECTOR,
    point, vector,
)

OUTPUT_LABEL   = "ShowerStartPosition"
DIRECTION_IN   = "ShowerDirection"

def run(
    vertex: Point | None,
    spacepoints: list[SpacePoint],
    holder: ShowerElementHolder,
    *,
    verbose: bool = False,
) -> int:
    """
    Calculate and store the shower start position.

    Parameters
    ----------
    vertex      : The PFP-associated interaction vertex (geo::Point_t equivalent).
                  Pass ``None`` to exercise the fall-back path.
    spacepoints : All space-points belonging to this shower PFParticle.
    holder      : The ShowerElementHolder shared across all tools.
    verbose     : If True, emit warning messages (mirrors fVerbose in C++).

    Returns
    -------
    0 on success, 1 on failure.
    """

    # use the PFP vertex
    if vertex is not None:
        err = UNSET_POINT.copy()
        holder.set_element(OUTPUT_LABEL, vertex.copy(), err)
        if verbose:
            print(f"[ShowerPFPVertexStartPosition] Start position set from PFP vertex: {vertex}")
        return 0

    # fall-back: order space-points along the already-computed PCA
    # direction and pick the most upstream point
    if holder.check_element(DIRECTION_IN):
        direction: Vector = holder.get_element(DIRECTION_IN)

        if len(spacepoints) == 0:
            if verbose:
                print("[ShowerPFPVertexStartPosition] No spacepoints available for fall-back.")
            return 0

        positions = numpy.array([sp.position for sp in spacepoints])  # (N, 3)

        # Shower centre
        centre = positions.mean(axis=0)

        # Project each point onto the shower axis (signed distance from centre)
        projections = (positions - centre) @ direction

        # Most negative projection = most upstream end
        idx_start = int(numpy.argmin(projections))
        start_pos = positions[idx_start].copy()

        err = UNSET_POINT.copy()
        holder.set_element(OUTPUT_LABEL, start_pos, err)

        if verbose:
            print(f"[ShowerPFPVertexStartPosition] Start position set from spacepoint ordering: "
                  f"{start_pos}")
        return 0

    # nothing we can do...
    if verbose:
        print("[ShowerPFPVertexStartPosition] Start position has not been set. "
              "No vertex and no direction available.")
        
    return 0