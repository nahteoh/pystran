"""
Define section dictionaries.

A section defines material properties, geometrical properties, such as the
second moment of area, and also orientation of the cross section profile.
"""

from math import pi
import numpy


def truss_section(name, E=0.0, A=0.0, rho=0.0, CTE=0.0):
    """
    Define truss section.

    - `E`= Young's modulus,
    - `A`= cross-sectional area,
    - `rho`= mass density,
    - `CTE`= coefficient of thermal expansion.
    """
    s = dict()
    s["name"] = name
    s["E"] = E
    s["rho"] = rho
    s["CTE"] = CTE
    s["A"] = A
    return s


def beam_2d_section(name, E=0.0, A=0.0, I=0.0, rho=0.0, CTE=0.0):
    """
    Define 2d beam section.

    - `E`= Young's modulus,
    - `A`= cross-sectional area,
    - `I`= central moment of inertia of the cross-section about the `y`
    coordinate axis (i.e. the axis perpendicular to the plane of the bending, `x-z`).
    - `rho`= mass density,
    - `CTE`= coefficient of thermal expansion.
    """
    s = dict()
    s["name"] = name
    s["E"] = E
    s["rho"] = rho
    s["CTE"] = CTE
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
    CTE=0.0,
):
    """
    Define 3d beam section.

    - `E`, `G`= Young's and shear modulus,
    - `rho`= mass density,
    - `A`= cross-sectional area,
    - `Ix`= central moment of inertia of the cross-section about the local `x`.
    - `Iy`, `Iz`= central moment of inertia of the cross-section about the local `y` and local `z`
    coordinate axis,
    - `J`= St Venant torsion constant.
    - `xz_vector`= vector that lies in the local `x` and `z` coordinate plane,
    - `CTE`= coefficient of thermal expansion.
    """
    s = dict()
    s["name"] = name
    s["E"] = E
    s["G"] = G
    s["rho"] = rho
    s["CTE"] = CTE
    s["A"] = A
    s["Ix"] = Ix
    s["Iy"] = Iy
    s["Iz"] = Iz
    s["J"] = J
    s["xz_vector"] = xz_vector
    return s


def circular_tube(innerradius, outerradius):
    """
    Returns the area, moments of inertia and torsion constant for a hollow circle (tube).
    """
    Rext = outerradius
    Rint = innerradius
    A = pi * (Rext ** 2 - Rint ** 2)
    Iy = pi / 4 * (Rext ** 4 - Rint ** 4)
    Iz = pi / 4 * (Rext ** 4 - Rint ** 4)
    Ix = Iy + Iz
    J = pi / 2 * (Rext ** 4 - Rint ** 4)
    return A, Ix, Iy, Iz, J


def i_beam(H, B, tf, tw):
    """
    Returns the area, moments of inertia and torsion constant for an I-beam.

    The axis parallel to the flanges is `y`, the axis parallel to the web is `z`.
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


def square_tube(H, B, th, tb):
    """
    Returns the area, moments of inertia and torsion constant for a square tube.

    The axis parallel to the `B` dimension is `y`, the axis parallel to the `H`
    dimension is `z`.
    """
    Bi, Hi = (B - 2 * tb), (H - 2 * th)
    A = B * H - Bi * Hi
    Iy = (B / 12) * H**3 - (Bi / 12) * Hi**3
    Iz = (B**3 / 12) * H - (Bi**3 / 12) * Hi
    Ix = Iy + Iz
    J = 2 * tb * th * Hi**2 * Bi**2 / (H * tb + B * th - tb**2 - th**2)
    return A, Ix, Iy, Iz, J


def rectangle(H, B):
    """
    Returns the area, moments of inertia and torsion constant for a solid
    rectangular section.

    The axis parallel to the `B` dimension is `y`, the axis parallel to the `H`
    dimension is `z`.
    """
    a = max(H, B)
    b = min(H, B)
    A = B * H
    Iy = (B / 12) * H**3
    Iz = (B**3 / 12) * H
    Ix = Iy + Iz
    rs = numpy.array([1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 10, 20, 40, 80, 200, 2000])
    coeff = numpy.array(
        [
            0.141,
            0.196,
            0.229,
            0.249,
            0.263,
            0.281,
            0.291,
            0.299,
            0.312,
            0.317,
            0.325,
            0.33,
            1 / 3,
            1 / 3,
        ]
    )
    c = numpy.interp(a / b, rs, coeff, coeff[0], coeff[-1])
    J = c * a * b**3
    return A, Ix, Iy, Iz, J
