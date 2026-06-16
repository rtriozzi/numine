"""
    showerreco.tools.shower_length
    ------------------------------
    Python equivalent of ShowerLengthPercentile.

    Logic:
    1. Project all spacepoints onto the shower axis (dot with unit direction
        from start position).  Sort by projection → pick the value at
        percentile P as the shower length.
    2. Project all spacepoints onto the perpendicular plane (transverse
        distance from axis).  Sort → pick value at percentile P as the width.
    3. Opening angle = atan(width / length).

    The percentile trick avoids outlier spacepoints dominating the length
    (typical value in icaruscode FHiCL: Percentile = 0.95).
"""

from __future__ import annotations

import numpy

from showerreco.geometry import ShowerElementHolder, SpacePoint


START_IN  = "ShowerStartPosition"
DIR_IN    = "ShowerDirection"
LENGTH_OUT = "ShowerLength"
ANGLE_OUT  = "ShowerOpeningAngle"


def run(
    spacepoints: list[SpacePoint],
    holder: ShowerElementHolder,
    *,
    percentile: float = 0.95,
    verbose: bool = False,
) -> int:
    """
        Calculate shower length and opening angle via longitudinal/transverse
        percentile cuts on the space-point distribution.

        Parameters
        ----------
        spacepoints : Space-points of the shower PFParticle.
        holder      : ShowerElementHolder (must have start position and direction).
        percentile  : Fraction of hits to include (mirrors fPercentile in C++).
        verbose     : Print diagnostics.

        Returns
        -------
        0 on success, 1 on failure.
    """

    if not holder.check_element(START_IN):
        if verbose:
            print("[ShowerLengthPercentile] Start position not set, returning.")
        return 1

    if not holder.check_element(DIR_IN):
        if verbose:
            print("[ShowerLengthPercentile] Direction not set, returning.")
        return 1

    if len(spacepoints) == 0:
        if verbose:
            print("[ShowerLengthPercentile] No spacepoints, returning.")
        return 1

    start = holder.get_element(START_IN)
    dirn  = holder.get_element(DIR_IN)
    u     = dirn / numpy.linalg.norm(dirn)  # ensure unit vector

    positions = numpy.array([sp.position for sp in spacepoints])    # (N, 3)
    delta     = positions - start   # (N, 3)

    # longitudinal projections (signed distance along shower axis)
    proj_long = delta @ u   # (N,)
    proj_long_sorted = numpy.sort(proj_long)

    length_idx    = int(numpy.floor(percentile * len(proj_long_sorted)))
    length_idx    = min(length_idx, len(proj_long_sorted) - 1)

    shower_length = float(proj_long_sorted[length_idx])
    shower_length_err = float(proj_long_sorted[-1] - shower_length) # max - percentile

    # transverse projection
    long_component = numpy.outer(proj_long, u)  # (N, 3)
    transverse     = delta - long_component # (N, 3)
    proj_perp      = numpy.linalg.norm(transverse, axis=1)  # (N,)
    proj_perp_sorted = numpy.sort(proj_perp)

    perp_idx   = int(numpy.floor(percentile * len(proj_perp_sorted)))
    perp_idx   = min(perp_idx, len(proj_perp_sorted) - 1)

    shower_width = float(proj_perp_sorted[perp_idx])

    # opening angle
    if shower_length > 0:
        shower_angle = float(numpy.arctan(shower_width / shower_length))
    else:
        shower_angle = 0.0

    shower_angle_err = -999.0   # not implemented (matches C++ TODO)

    # store
    holder.set_element(LENGTH_OUT, shower_length, shower_length_err)
    holder.set_element(ANGLE_OUT,  shower_angle,  shower_angle_err)

    if verbose:
        print(f"[ShowerLengthPercentile] Length : {shower_length:.2f} cm  "
              f"(err +{shower_length_err:.2f} cm)")
        print(f"[ShowerLengthPercentile] Width  : {shower_width:.2f} cm")
        print(f"[ShowerLengthPercentile] Angle  : {numpy.degrees(shower_angle):.2f} deg")

    return 0