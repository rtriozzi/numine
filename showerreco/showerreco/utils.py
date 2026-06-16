# showerreco/utils.py

import numpy as numpy

import matplotlib 
import matplotlib.pyplot as plt

"""
    Cone overlay for the shower 3D event display.
    Usage:
        add_shower_cone(ax, START, DIRECTION, LENGTH, ANGLE, n_rings=30, n_phi=60)
"""

def add_3D_shower_cone(
    ax,
    start:     "array-like",   # [x, y, z] of shower start
    direction: "array-like",   # unit vector [dx, dy, dz] (PCA direction)
    length:    float,          # shower length in cm
    angle:     float,          # half-opening angle in radians
    *,
    n_rings:   int   = 40,     # longitudinal resolution
    n_phi:     int   = 80,     # azimuthal resolution
    color:     str   = "red",
    alpha:     float = 0.10,   # surface transparency
    wire_alpha: float = 0.30,  # wireframe transparency
    n_wires:   int   = 8,      # number of azimuthal wireframe lines
    label:     str   = "shower cone",
    ax_order:  tuple = (0, 2, 1),  # (x, y, z) → scatter order, default: X Z Y
):
    """
        Draw a semi-transparent cone on an existing Axes3D.

        The cone runs from `start` to `start + length * direction`.
        Its radius at longitudinal position s is  r(s) = s * tan(angle).

        Parameters
        ----------
        ax_order : tuple of 3 ints
            Index permutation applied when calling ax.plot_surface / ax.plot.
            Your existing scatter uses (X, Z, Y), so the default is (0, 2, 1).
            Change to (0, 1, 2) if your axes are in the natural X Y Z order.
    """
    start     = numpy.asarray(start,     dtype=float)
    direction = numpy.asarray(direction, dtype=float)
    direction = direction / numpy.linalg.norm(direction)

    # build an orthonormal basis {u, v} perpendicular to direction
    ref = numpy.array([0., 0., 1.]) if abs(direction[2]) < 0.9 else numpy.array([1., 0., 0.])
    u   = numpy.cross(direction, ref);  u /= numpy.linalg.norm(u)
    v   = numpy.cross(direction, u)

    # parametric surface:  P(s, phi) = start + s*direction + r(s)*(cos(phi)*u + sin(phi)*v)
    s_vals   = numpy.linspace(0., length, n_rings)
    phi_vals = numpy.linspace(0., 2 * numpy.pi, n_phi)
    S, PHI   = numpy.meshgrid(s_vals, phi_vals)

    R = S * numpy.tan(angle)   # radius grows linearly

    # 3D coordinates of the cone surface
    P = (
        start[None, None, :]
        + S[:, :, None]   * direction[None, None, :]
        + R[:, :, None]   * (numpy.cos(PHI)[:, :, None] * u[None, None, :]
                           +  numpy.sin(PHI)[:, :, None] * v[None, None, :])
    )  # shape (n_phi, n_rings, 3)

    i0, i1, i2 = ax_order   # e.g. (0, 2, 1) → X, Z, Y

    # filled surface
    ax.plot_surface(
        P[:, :, i0], P[:, :, i1], P[:, :, i2],
        alpha     = alpha,
        color     = color,
        linewidth = 0,
        antialiased = True,
    )

    # rim circle at the end
    phi_rim = numpy.linspace(0., 2 * numpy.pi, n_phi)
    r_rim   = length * numpy.tan(angle)
    rim     = (
        start[None, :]
        + length * direction[None, :]
        + r_rim * (numpy.cos(phi_rim)[:, None] * u[None, :]
                 + numpy.sin(phi_rim)[:, None] * v[None, :])
    )
    ax.plot(
        rim[:, i0], rim[:, i1], rim[:, i2],
        color=color, linewidth=1.0, alpha=wire_alpha,
    )

    # # azimuthal wireframe lines from apex to rim
    # for phi_w in numpy.linspace(0., 2 * numpy.pi, n_wires, endpoint=False):
    #     rim_pt = (
    #         start
    #         + length * direction
    #         + r_rim * (numpy.cos(phi_w) * u + numpy.sin(phi_w) * v)
    #     )
    #     ax.plot(
    #         [start[i0], rim_pt[i0]],
    #         [start[i1], rim_pt[i1]],
    #         [start[i2], rim_pt[i2]],
    #         color=color, linewidth=0.8, alpha=wire_alpha,
    #     )

    # invisible proxy for the legend
    proxy = plt.matplotlib.lines.Line2D(
        [], [], color=color, linewidth=2, alpha=wire_alpha, label=label,
    )
    ax.add_artist(proxy)
    return proxy

def rotate_to_vector(
    v, 
    target
):
    target = target / numpy.linalg.norm(target)
    v = v / numpy.linalg.norm(v)

    axis = numpy.cross(v, target)
    s = numpy.linalg.norm(axis)
    c = numpy.dot(v, target)

    if s < 1e-8:
        return numpy.eye(3)

    axis = axis / s

    K = numpy.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])

    R = numpy.eye(3) + K + K @ K * ((1 - c) / (s**2))
    
    return R