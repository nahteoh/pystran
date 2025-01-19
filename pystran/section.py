"""
Define section dictionaries.

A section defines material properties, geometrical properties, such as the
second moment of area, and also orientation of the cross section profile.
"""

from math import pi
import numpy


def truss_section(name, E=0.0, A=0.0, rho=0.0):
    """
    Define truss section.
    """
    s = dict()
    s["name"] = name
    s["E"] = E
    s["rho"] = rho
    s["A"] = A
    return s


def beam_2d_section(name, E=0.0, A=0.0, I=0.0, rho=0.0):
    """
    Define 2d beam section.

    - `E`= Young's modulus,
    - `A`= cross-sectional area,
    - `I`= central moment of inertia of the cross-section about the x3
    coordinate axis (i.e. the axis perpendicular to the plane of the beam).
    """
    s = dict()
    s["name"] = name
    s["E"] = E
    s["rho"] = rho
    s["A"] = A
    s["I"] = I
    return s


def beam_3d_section(
    name,
    E=0.0,
    G=0.0,
    A=0.0,
    Ix=0.0,
    Iy=0.0,
    Iz=0.0,
    J=0.0,
    rho=0.0,
    xz_vector=(0, 0, 1),
):
    """
    Define 3d beam section.

    - `E`, `G`= Young's and shear modulus,
    - `A`= cross-sectional area,
    - `Ix`= central moment of inertia of the cross-section about the local x.
    - `Iy`, `Iz`= central moment of inertia of the cross-section about the local y and local z
    coordinate axis,
    - `J`= St Venant torsion constant.
    - `xz_vector`= vector that lies in the local x and z coordinate plane.
    """
    s = dict()
    s["name"] = name
    s["E"] = E
    s["G"] = G
    s["rho"] = rho
    s["A"] = A
    s["Ix"] = Ix
    s["Iy"] = Iy
    s["Iz"] = Iz
    s["J"] = J
    s["xz_vector"] = xz_vector
    return s


def close_points(points):
    """
    If the input points do not form a closed polygon, closes the polygon
    and returns the result.

    Parameters
    ----------
    points : array
        An array of (x, y) coordinates of shape (N, 2).
    """
    p = numpy.asarray(points)
    if (p[0] == p[-1]).all():
        return p
    return numpy.append(p, [p[0]], axis=0)


def hollow_circle(innerradius, outerradius):
    """
    Returns the area, moments of inertia and torsion constant for a hollow circle (tube).
    """
    Rext = outerradius
    Rint = innerradius
    A = pi * (Rext ^ 2 - Rint ^ 2)
    Iy = pi / 4 * (Rext ^ 4 - Rint ^ 4)
    Iz = pi / 4 * (Rext ^ 4 - Rint ^ 4)
    Ix = Iy + Iz
    J = pi / 2 * (Rext ^ 4 - Rint ^ 4)
    return A, Ix, Iy, Iz, J


def i_beam(H, B, tf, tw):
    """
    Returns the area, moments of inertia and torsion constant for an I-beam.
    """
    A = B * H - (B - tw) * (H - 2 * tf)
    Iy = (B / 12) * H**3 - ((B - tw) / 12) * (H - 2 * tf) ** 3
    Iz = (
        H * B**3 / 12
        - 2 * ((B - tw) / 2) ** 3 * (H - 2 * tf) / 12
        - 2 * ((B - tw) / 2) * (H - 2 * tf) * ((B - tw) / 4 + tw / 2) ** 2
    )
    Ix = Iy + Iz
    J = (2 * B * tf**3 + (H - 2 * tf) * tw**3) / 3
    return A, Ix, Iy, Iz, J
