"""
    showerreco.tools.pca_direction
    ------------------------------
    Python equivalent of ShowerPCADirection.

    Original C++ logic (ShowerPCADirection.cc):
    1. Retrieve the space-points for the shower PFParticle.
    2. Optionally charge-weight them (correcting for electron lifetime).
    3. Build the 3×3 covariance matrix and diagonalise it (Eigen).
    4. The primary eigenvector (largest eigenvalue) is the raw PCA direction.
    5. Orient the vector so that it points *away* from the start position
        (dot product with start→centre vector must be positive).
    6. Fall back to RMS-gradient orientation if no start position is available.

    Electron-lifetime correction
    -----------------------------
    In the C++ code the charge of each space-point is corrected before weighting:

        Q_corr = Q * exp( (sampling_rate * t) / (tau_e * 1e3) )

    where tau_e is in ms and sampling_rate * t gives the drift time in µs.
    We keep the same formula; caller can set electron_lifetime_ms=None to skip.
    """

from __future__ import annotations

import numpy

from showerreco.geometry import (
    Point, Vector, SpacePoint, PCAxis,
    ShowerElementHolder,
    UNSET_POINT, UNSET_VECTOR,
)

# element-holder keys
START_POS_IN    = "ShowerStartPosition"
DIRECTION_OUT   = "ShowerDirection"
CENTRE_OUT      = "ShowerCentre"
PCA_OUT         = "ShowerPCA"

# minimum number of space-points required to run PCA (same as C++)
MIN_SPACEPOINTS = 3

def run(
    spacepoints: list[SpacePoint],
    holder: ShowerElementHolder,
    *,
    charge_weighted: bool = True,
    use_start_position: bool = True,
    n_segments: int = 6,
    electron_lifetime_ms: float | None = 3.0,   # ms; None → skip correction
    sampling_rate_us: float = 0.5,               # µs per tick
    verbose: bool = False,
) -> int:
    """
        Compute the shower direction via PCA on the shower's space-points.

        Parameters
        ----------
        spacepoints          : Space-points belonging to this shower PFParticle.
        holder               : Shared ShowerElementHolder.
        charge_weighted      : Apply charge-weighting to the covariance matrix.
        use_start_position   : Orient PCA vector using start→centre direction.
        n_segments           : Number of segments for RMS-gradient orientation
                            (only used when use_start_position is False).
        electron_lifetime_ms : Electron lifetime in ms for charge correction.
        sampling_rate_us     : Detector sampling rate in µs/tick.
        verbose              : Print diagnostic messages.

        Returns
        -------
        0 on success, 1 on failure.
    """

    if len(spacepoints) < MIN_SPACEPOINTS:
        if verbose:
            print(f"[ShowerPCADirection] Only {len(spacepoints)} spacepoints — skipping.")
        return 1

    # compute the PCA
    pca, centre = _calculate_pca(
        spacepoints,
        charge_weighted=charge_weighted,
        electron_lifetime_ms=electron_lifetime_ms,
        sampling_rate_us=sampling_rate_us,
    )

    pca_direction: Vector = pca.eigenvectors[0].copy()   # primary eigenvector

    # store centre and full PCA object
    holder.set_element(CENTRE_OUT, centre.copy(), UNSET_POINT.copy())
    holder.set_element(PCA_OUT,    pca)

    # orient the direction vector
    if use_start_position:
        if not holder.check_element(START_POS_IN):
            if verbose:
                print("[ShowerPCADirection] use_start_position=True but start position not set.")
            return 1

        start: Point = holder.get_element(START_POS_IN)
        general_dir  = centre - start
        norm = numpy.linalg.norm(general_dir)
        if norm > 0:
            general_dir = general_dir / norm

        if numpy.dot(pca_direction, general_dir) < 0:
            pca_direction = -pca_direction

    else:
        # RMS-gradient method: gradient > 0 means the vector already points
        # from narrow end to wide end (downstream), keep it; otherwise flip.
        gradient = _rms_gradient(spacepoints, centre, pca_direction, n_segments)
        if verbose:
            print(f"[ShowerPCADirection] RMS gradient = {gradient:.4f}")
        if gradient < -numpy.finfo(float).eps:
            pca_direction = -pca_direction

    # store result
    holder.set_element(DIRECTION_OUT, pca_direction, UNSET_POINT.copy())

    if verbose:
        print(f"[ShowerPCADirection] Direction: {pca_direction}")
        print(f"[ShowerPCADirection] Centre:    {centre}")
        evals = pca.eigenvalues
        print(f"[ShowerPCADirection] Eigenvalues: {evals[0]:.3f}, {evals[1]:.3f}, {evals[2]:.3f}")

    return 0

# PCA
def _calculate_pca(
    spacepoints: list[SpacePoint],
    *,
    charge_weighted: bool,
    electron_lifetime_ms: float | None,
    sampling_rate_us: float,
) -> tuple[PCAxis, Point]:
    """Build the covariance matrix and return the PCAxis and shower centre."""

    positions = numpy.array([sp.position for sp in spacepoints], dtype=numpy.float64)  # (N, 3)

    # weights
    if charge_weighted:
        charges = numpy.array([sp.charge for sp in spacepoints], dtype=numpy.float64)
        times   = numpy.array([sp.time   for sp in spacepoints], dtype=numpy.float64)

        # electron-lifetime correction
        if electron_lifetime_ms is not None and electron_lifetime_ms > 0:
            tau_us = electron_lifetime_ms * 1e3          # convert ms → µs
            drift_us = sampling_rate_us * times
            charges = charges * numpy.exp(drift_us / tau_us)

        # guard against non-positive charges
        charges = numpy.clip(charges, 0.0, None)
        total_charge = charges.sum()
        if total_charge <= 0:
            charge_weighted = False         # fall back gracefully

    if charge_weighted:
        # charge-weighted centroid
        centre = (positions * charges[:, None]).sum(axis=0) / total_charge
        # sqrt-charge weights (same as C++)
        weights = numpy.sqrt(charges / total_charge)
    else:
        centre  = positions.mean(axis=0)
        weights = numpy.ones(len(spacepoints))

    # covariance matrix
    delta = positions - centre                         # (N, 3)
    w_delta = delta * weights[:, None]                 # (N, 3) weighted
    cov = (w_delta.T @ delta) / weights.sum()          # (3, 3)

    # diagonalise 
    eigenvalues, eigenvectors = numpy.linalg.eigh(cov)   # ascending order

    # reverse to descending (match Eigen's SelfAdjointEigenSolver output order)
    eigenvalues  = eigenvalues[::-1]
    eigenvectors = eigenvectors[:, ::-1].T             # rows = eigenvectors

    pca = PCAxis(
        ok=True,
        n_hits=len(spacepoints),
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        ave_position=centre,
    )
    return pca, centre

# RMS-gradient orientation
def _rms_gradient(
    spacepoints: list[SpacePoint],
    centre: Point,
    direction: Vector,
    n_segments: int,
) -> float:
    """
        Split the shower into n_segments along `direction`, compute the transverse
        RMS in each, and return the slope of a linear fit (RMS vs segment index).

        A positive slope means the shower widens going forward → direction is
        already correct.  A negative slope means flip it.
    """
    if n_segments < 2:
        return 0.0

    positions = numpy.array([sp.position for sp in spacepoints], dtype=numpy.float64)
    delta     = positions - centre

    # longitudinal projections
    projections = delta @ direction

    lo, hi = projections.min(), projections.max()
    if hi <= lo:
        return 0.0

    # transverse displacement magnitude
    long_component  = numpy.outer(projections, direction)
    transverse      = delta - long_component
    trans_mag       = numpy.linalg.norm(transverse, axis=1)

    # bin edges
    edges = numpy.linspace(lo, hi, n_segments + 1)
    seg_centres = 0.5 * (edges[:-1] + edges[1:])

    rms_values = []
    for i in range(n_segments):
        mask = (projections >= edges[i]) & (projections < edges[i + 1])
        if mask.sum() < 2:
            rms_values.append(0.0)
        else:
            rms_values.append(float(numpy.std(trans_mag[mask])))

    # linear fit: gradient of RMS vs segment index
    seg_indices = numpy.arange(n_segments, dtype=float)
    if numpy.all(numpy.array(rms_values) == 0):
        return 0.0

    coeffs = numpy.polyfit(seg_centres, rms_values, 1)
    return float(coeffs[0])   # slope