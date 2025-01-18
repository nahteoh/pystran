"""
Define beam mechanical quantities.
"""

from math import sqrt
from numpy import array, dot, outer, concatenate, zeros
from numpy.linalg import norm, cross
from pystran import geometry
from pystran import assemble
from pystran import truss


def beam_2d_shape_fun(xi):
    """
    Compute the beam shape functions for deflection in the x-z plane (i.e. in 2d).
    """
    return array(
        [
            (2 - 3 * xi + xi**3) / 4,
            (-1 + xi + xi**2 - xi**3) / 4,
            (2 + 3 * xi - xi**3) / 4,
            (+1 + xi - xi**2 - xi**3) / 4,
        ]
    )


def beam_2d_shape_fun_xi(xi):
    """
    Compute the first derivative of the beam shape functions for deflection in the x-z plane (i.e. in 2d).

    The quantity computed is
    ```math
    \frac{d^1 N(\\xi)}{d \\xi^1}
    ```
    """
    return array(
        [
            (-3 + 3 * xi**2) / 4,
            (+1 + 2 * xi - 3 * xi**2) / 4,
            (3 - 3 * xi**2) / 4,
            (+1 - 2 * xi - 3 * xi**2) / 4,
        ]
    )


def beam_2d_shape_fun_xi2(xi):
    """
    Compute the second derivative of the beam shape functions for deflection in the x-z plane (i.e. in 2d).

    The quantity computed is
    ```math
    \frac{d^2 N(\\xi)}{d \\xi^2}
    ```
    """
    return array([(6 * xi) / 4, (2 - 6 * xi) / 4, (-6 * xi) / 4, (-2 - 6 * xi) / 4])


def beam_2d_shape_fun_xi3(xi):
    """
    Compute the third derivative of the beam shape functions for deflection in the x-z plane (i.e. in 2d).

    The quantity computed is
    ```math
    \frac{d^3 N(\\xi)}{d \\xi^3}
    ```
    """
    return array([(6) / 4, (-6) / 4, (-6) / 4, (-6) / 4])


def beam_3d_xz_shape_fun(xi):
    """
    Compute the beam shape functions for deflection in the x-z plane.
    """
    return beam_2d_shape_fun(xi)


def beam_3d_xz_shape_fun_xi2(xi):
    """
    Compute the second derivative of the beam shape functions for deflection in the x-z plane.
    """
    return beam_2d_shape_fun_xi2(xi)


def beam_3d_xy_shape_fun(xi):
    """
    Compute the beam shape functions for deflection in the x-y plane.

    The quantity computed is
    ```math
    \frac{d^2 N(\\xi)}{d \\xi^2}
    ```

    The signs of the shape functions that go with the rotations (i.e. the second and fourth) need to be reversed.
    """
    N = beam_2d_shape_fun(xi)
    N[1] *= -1.0
    N[3] *= -1.0
    return N


def beam_3d_xy_shape_fun_xi2(xi):
    """
    Compute the second derivative of the beam shape functions for deflection in the x-y plane.

    The quantity computed is
    ```math
    \frac{d^2 N(\\xi)}{d \\xi^2}
    ```

    The signs of the shape functions that go with the rotations (i.e. the second and fourth) need to be reversed.
    """
    d2Ndxi2 = beam_2d_shape_fun_xi2(xi)
    d2Ndxi2[1] *= -1.0
    d2Ndxi2[3] *= -1.0
    return d2Ndxi2


def beam_2d_member_geometry(i, j):
    """
    Compute 2d beam geometry.

    The deformation of the beam is considered in the x-z plane. `e_x` is the
    direction vector along the axis of the beam. `e_z` is the direction vector
    perpendicular to the axis of the beam. These two vectors form a right-handed
    coordinate system, considering `e_y` points out of the whiteboard
    (consistent with the sign convention in the book).
    """
    e_x = geometry.delt(i["coordinates"], j["coordinates"])
    h = geometry.len(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h
    # The orientation here reflects the sign convention in the book.
    # The deflection is measured positive downwards, while the x coordinate is measured left to right.
    # So in two dimensions e_x and e_z form a left-handed coordinate system.
    e_z = array([e_x[1], -e_x[0]])
    return e_x, e_z, h


def beam_3d_member_geometry(i, j, xz_vector):
    """
    Compute 3d beam geometry.
    """
    e_x = geometry.delt(i["coordinates"], j["coordinates"])
    h = geometry.len(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h
    # The orientation here reflects the sign convention in the book.
    # The deflection is measured positive downwards, while the x coordinate is measured left to right.
    # So in two dimensions e_x and e_z form a left-handed coordinate system.
    if abs(dot(e_x, xz_vector)) > 0.99 * norm(xz_vector):
        raise ZeroDivisionError("xz_vector must not be parallel to the beam axis")
    e_y = cross(xz_vector, e_x)
    e_y = e_y / norm(e_y)
    e_z = cross(e_x, e_y)
    return e_x, e_y, e_z, h


def beam_2d_bending_stiffness(e_x, e_z, h, E, I):
    """
    Compute 2d beam stiffness matrix.
    """
    xiG = [-1 / sqrt(3), 1 / sqrt(3)]
    WG = [1, 1]
    K = zeros((6, 6))
    for q in range(2):
        B = beam_2d_curvature_displacement(e_x, e_z, h, xiG[q])
        K += E * I * outer(B.T, B) * WG[q] * (h / 2)
    return K


def beam_3d_xz_curvature_displacement(e_x, e_y, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix in the local x-z plane.

    The quantity computed is
    ```math
    \\frac{d^2 N(x)}{d x^2} =  \\frac{d^2 N(\\xi)}{d \\xi^2} (2/h)^2
    ```
    """
    d2Ndxi2 = beam_2d_shape_fun_xi2(xi)
    B = zeros((1, 12))
    B[0, 0:3] = d2Ndxi2[0] * (2 / h) ** 2 * e_z
    B[0, 3:6] = (h / 2) * d2Ndxi2[1] * (2 / h) ** 2 * e_y
    B[0, 6:9] = d2Ndxi2[2] * (2 / h) ** 2 * e_z
    B[0, 9:12] = (h / 2) * d2Ndxi2[3] * (2 / h) ** 2 * e_y
    return B


def beam_3d_xy_curvature_displacement(e_x, e_y, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix in the local x-y plane.

    The quantity computed is
    ```math
    \\frac{d^2 N(x)}{d x^2} =  \\frac{d^2 N(\\xi)}{d \\xi^2} (2/h)^2
    ```
    """
    B = zeros((1, 12))
    B[0, 0:3] = 6 * xi / h**2 * e_y
    B[0, 3:6] = -(1 - 3 * xi) / h * e_z
    B[0, 6:9] = -6 * xi / h**2 * e_y
    B[0, 9:12] = +(3 * xi + 1) / h * e_z
    return B


def beam_3d_bending_stiffness(e_x, e_y, e_z, h, E, Iy, Iz):
    """
    Compute 2d beam stiffness matrix.
    """
    xiG = [-1 / sqrt(3), 1 / sqrt(3)]
    WG = [1, 1]
    Kxy = zeros((12, 12))
    for q in range(2):
        B = beam_3d_xy_curvature_displacement(e_x, e_y, e_z, h, xiG[q])
        Kxy += E * Iz * outer(B.T, B) * WG[q] * (h / 2)
    Kxz = zeros((12, 12))
    for q in range(2):
        B = beam_3d_xz_curvature_displacement(e_x, e_y, e_z, h, xiG[q])
        Kxz += E * Iy * outer(B.T, B) * WG[q] * (h / 2)
    return Kxy, Kxz


def beam_2d_curvature_displacement(e_x, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix.
    """
    d2Ndxi2 = beam_2d_shape_fun_xi2(xi)
    B = zeros((1, 6))
    B[0, 0:2] = d2Ndxi2[0] * (2 / h) ** 2 * e_z
    B[0, 2] = (h / 2) * d2Ndxi2[1] * (2 / h) ** 2
    B[0, 3:5] = d2Ndxi2[2] * (2 / h) ** 2 * e_z
    B[0, 5] = (h / 2) * d2Ndxi2[3] * (2 / h) ** 2
    return B


def beam_2d_3rd_deriv_displacement(e_x, e_z, h, xi):
    """
    Compute beam third derivative-displacement matrix.
    """
    d2Ndxi3 = beam_2d_shape_fun_xi3(xi)
    B = zeros((1, 6))
    B[0, 0:2] = d2Ndxi3[0] * (2 / h) ** 3 * e_z
    B[0, 2] = (h / 2) * d2Ndxi3[1] * (2 / h) ** 3
    B[0, 3:5] = d2Ndxi3[2] * (2 / h) ** 3 * e_z
    B[0, 5] = (h / 2) * d2Ndxi3[3] * (2 / h) ** 3
    return B


def beam_2d_moment(member, i, j, xi):
    """
    Compute 2d beam moment based on the displacements stored at the joints.
    The moment is computed at the parametric location `xi` along the beam.
    """
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    sect = member["section"]
    E, I = sect["E"], sect["I"]
    ui, uj = i["displacements"], j["displacements"]
    u = concatenate([ui, uj])
    B = beam_2d_curvature_displacement(e_x, e_z, h, xi)
    return -E * I * dot(B, u)


def beam_3d_moment(member, i, j, axis, xi):
    """
    Compute 3d beam moment based on the displacements stored at the joints.
    The moment is computed at the parametric location `xi` along the beam.
    """
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    E, Iy, Iz = sect["E"], sect["Iy"], sect["Iz"]
    ui, uj = i["displacements"], j["displacements"]
    u = concatenate([ui, uj])
    if axis == "y":
        B = beam_3d_xz_curvature_displacement(e_x, e_y, e_z, h, xi)
        M = -E * Iy * dot(B, u)
    else:
        B = beam_3d_xy_curvature_displacement(e_x, e_y, e_z, h, xi)
        M = -E * Iz * dot(B, u)
    return M


def beam_2d_shear_force(member, i, j, xi):
    """
    Compute 2d beam shear force based on the displacements stored at the joints.
    The shear force is computed at the parametric location `xi` along the beam.
    """
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    sect = member["section"]
    E, I = sect["E"], sect["I"]
    ui, uj = i["displacements"], j["displacements"]
    u = concatenate([ui, uj])
    B = beam_2d_3rd_deriv_displacement(e_x, e_z, h, xi)
    return -E * I * dot(B, u)


def _stiffness_truss(member, i, j):
    e_x, L = truss.truss_member_geometry(i, j)
    sect = member["section"]
    E, A = sect["E"], sect["A"]
    return truss.stiffness(e_x, L, E, A)


def torsion_stiffness(e_x, h, G, J):
    """
    Compute torsion stiffness matrix.
    """
    B = torsion_displacement(e_x, h)
    return G * J * outer(B.T, B) * h


def torsion_displacement(e_x, h):
    """
    Compute torsion-displacement matrix.
    """
    B = zeros((1, 6))
    B[0, 0:3] = -e_x / h
    B[0, 3:6] = e_x / h
    return B


def _stiffness_torsion(member, i, j):
    e_x, h = truss.truss_member_geometry(i, j)
    sect = member["section"]
    G, J = sect["G"], sect["J"]
    return torsion_stiffness(e_x, h, G, J)


def assemble_stiffness(Kg, member, i, j):
    """
    Assemble beam stiffness matrix.
    """
    # Add stiffness in bending.
    beam_is_2d = len(i["coordinates"]) == len(j["coordinates"]) == 2
    dof = concatenate([i["dof"], j["dof"]])
    if beam_is_2d:
        e_x, e_z, h = beam_2d_member_geometry(i, j)
        sect = member["section"]
        E, I = sect["E"], sect["I"]
        k = beam_2d_bending_stiffness(e_x, e_z, h, E, I)
        Kg = assemble.assemble(Kg, dof, k)
    else:
        sect = member["section"]
        e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
        E, Iy, Iz = sect["E"], sect["Iy"], sect["Iz"]
        kxy, kxz = beam_3d_bending_stiffness(e_x, e_y, e_z, h, E, Iy, Iz)
        Kg = assemble.assemble(Kg, dof, kxy)
        Kg = assemble.assemble(Kg, dof, kxz)
    # Add stiffness in the axial direction.
    k = _stiffness_truss(member, i, j)
    if beam_is_2d:
        dof = concatenate([i["dof"][0:2], j["dof"][0:2]])
    else:
        dof = concatenate([i["dof"][0:3], j["dof"][0:3]])
    Kg = assemble.assemble(Kg, dof, k)
    if not beam_is_2d:
        # Add stiffness in torsion.
        k = _stiffness_torsion(member, i, j)
        dof = concatenate([i["dof"][3:6], j["dof"][3:6]])
        Kg = assemble.assemble(Kg, dof, k)
    return Kg
