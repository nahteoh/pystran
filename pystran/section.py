"""
Define section dictionaries.

A section defines material properties, geometrical properties, such as the
second moment of area, and also orientation of the cross section profile.
"""


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
    xz_vector=[0, 0, 1],
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
