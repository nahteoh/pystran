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
    Compute the beam shape functions for deflection in the `x-z` plane (i.e. in 2d).

    An array of shape function values is returned (i.e. $[N_1(\\xi), ..., N_4(\\xi)]$).
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
    Compute the first derivative of the beam shape functions for deflection in
    the `x-z` plane (i.e. in 2d).

    An array of first derivatives of shape functions is returned (i.e.
    $[dN_1(\\xi)/d\\xi, ..., dN_4(\\xi)/d\\xi]$).
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
    Compute the second derivative of the beam shape functions for deflection in
    the `x-z` plane (i.e. in 2d).

    An array of second derivatives of shape functions is returned (i.e.
    $[d^2N_1(\\xi)/d\\xi^2, ..., d^2N_4(\\xi)/d\\xi^2]$).
    """
    return array([(6 * xi) / 4, (2 - 6 * xi) / 4, (-6 * xi) / 4, (-2 - 6 * xi) / 4])


def beam_2d_shape_fun_xi3(xi):
    """
    Compute the third derivative of the beam shape functions for deflection in
    the `x-z` plane (i.e. in 2d).

    An array of third derivatives of shape functions is returned (i.e.
    $[d^3N_1(\\xi)/d\\xi^3, ..., d^3N_4(\\xi)/d\\xi^3]$).
    """
    return array([(6) / 4, (-6) / 4, (-6) / 4, (-6) / 4])


def beam_3d_xz_shape_fun(xi):
    """
    Compute the beam shape functions for deflection in the `x-z` plane.

    An array of shape function values is returned (i.e. $[N_1(\\xi), ..., N_4(\\xi)]$).
    """
    return beam_2d_shape_fun(xi)


def beam_3d_xz_shape_fun_xi2(xi):
    """
    Compute the second derivative of the beam shape functions for deflection in
    the `x-z` plane.

    An array of second derivatives of shape functions is returned (i.e.
    $[d^2N_1(\\xi)/d\\xi^2, ..., d^2N_4(\\xi)/d\\xi^2]$).
    """
    return beam_2d_shape_fun_xi2(xi)


def beam_3d_xy_shape_fun(xi):
    """
    Compute the beam shape functions for deflection in the x-y plane.

    The signs of the shape functions that go with the rotations (i.e. the second
    and fourth) need to be reversed: An array of second derivatives of shape
    functions is returned (i.e. $[N_1(\\xi), -N_2(\\xi),
    N_3(\\xi), -N_4(\\xi)]$).
    """
    N = beam_2d_shape_fun(xi)
    N[1] *= -1.0
    N[3] *= -1.0
    return N


def beam_3d_xy_shape_fun_xi2(xi):
    """
    Compute the second derivative of the beam shape functions for deflection in
    the `x-y` plane.

    The signs of the shape functions that go with the rotations (i.e. the second
    and fourth) need to be reversed: An array of second derivatives of shape
    functions is returned (i.e. $[d^2N_1(\\xi)/d\\xi^2, -d^2N_2(\\xi)/d\\xi^2,
    d^2N_3(\\xi)/d\\xi^2, -d^2N_4(\\xi)/d\\xi^2]$).
    """
    d2Ndxi2 = beam_2d_shape_fun_xi2(xi)
    d2Ndxi2[1] *= -1.0
    d2Ndxi2[3] *= -1.0
    return d2Ndxi2


def beam_3d_xz_shape_fun_xi3(xi):
    """
    Compute the third derivative of the beam shape functions with respect to
    `xi` for deflection in the `x-z` plane.

    An array of third derivatives of shape functions is returned (i.e.
    $[d^3N_1(\\xi)/d\\xi^3, ..., d^3N_4(\\xi)/d\\xi^3]$).
    """
    return beam_2d_shape_fun_xi3(xi)


def beam_3d_xy_shape_fun_xi3(xi):
    """
    Compute the third derivative of the beam shape functions with respect to
    `xi` for deflection in the `x-y` plane.

    The signs of the shape functions that go with the rotations (i.e. the second
    and fourth) need to be reversed: An array of third derivatives of shape
    functions is returned (i.e. $[d^3N_1(\\xi)/d\\xi^3, -d^3N_2(\\xi)/d\\xi^3,
    d^3N_3(\\xi)/d\\xi^3, -d^3N_4(\\xi)/d\\xi^3]$).
    """
    d3Ndxi3 = beam_2d_shape_fun_xi3(xi)
    d3Ndxi3[1] *= -1.0
    d3Ndxi3[3] *= -1.0
    return d3Ndxi3


def beam_2d_member_geometry(i, j):
    """
    Compute 2d beam geometry.

    A local coordinate system is attached to the beam such that the `x` axis is
    along the beam axis. The deformation of the beam is considered in the x-z
    plane.

    Vector `e_x` is the direction vector along the axis of the beam. `e_z` is
    the direction vector perpendicular to the axis of the beam. These two
    vectors form a left-handed coordinate system (consistent with the sign
    convention in the book): The deflection `w` is measured positive downwards,
    while the `x` coordinate is measured left to right. So in two dimensions
    `e_x` and `e_z` form a left-handed coordinate system. In reality, the
    complete coordinate system is right-handed, as the not-used basis vector is
    `e_y`, which points out of the plane of the screen (page).
    """
    e_x = geometry.delt(i["coordinates"], j["coordinates"])
    h = geometry.vlen(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h
    e_z = array([e_x[1], -e_x[0]])
    return e_x, e_z, h


def beam_3d_member_geometry(i, j, xz_vector):
    """
    Compute 3d beam geometry.

    A local coordinate system is attached to the beam such that the `x` axis is
    along the beam axis. The deformation of the beam is considered in the `x-y`
    and `x-z` plane.

    Vector `e_x` is the direction vector along the axis of the beam. `e_z` is
    the direction vector perpendicular to the axis of the beam. These two
    vectors form a right-handed coordinate system, completed by `e_y`.

    The plane `x-z` is defined by the vector `xz_vector` and the beam axis (i.e.
    `e_x`). Therefore, the vector `xz_vector` must not be parallel to the beam
    axis.

    - `i` and `j` = the two joints of the beam.
    - `xz_vector` = the vector that defines the `x-z` plane of the beam-local
      coordinate system. It does not need to be of unit length, but it must not
      be parallel to the beam axis.
    """
    e_x = geometry.delt(i["coordinates"], j["coordinates"])
    h = geometry.vlen(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h  # normalize the unit length
    # The orientation here reflects the sign convention in the book. The
    # deflection is measured positive downwards, while the x coordinate is
    # measured left to right. So in two dimensions e_x and e_z form a
    # left-handed coordinate system.
    if abs(dot(e_x, xz_vector)) > 0.99 * norm(xz_vector):
        raise ZeroDivisionError("xz_vector must not be parallel to the beam axis")
    e_y = cross(xz_vector, e_x)
    e_y = e_y / norm(e_y)
    e_z = cross(e_x, e_y)
    return e_x, e_y, e_z, h


def beam_2d_bending_stiffness(e_z, h, E, I):
    """
    Compute 2d beam stiffness matrix.
    """
    xiG = [-1 / sqrt(3), 1 / sqrt(3)]
    WG = [1, 1]
    K = zeros((6, 6))
    for q in range(2):
        B = beam_2d_curv_displ_matrix(e_z, h, xiG[q])
        K += E * I * outer(B.T, B) * WG[q] * (h / 2)
    return K


def beam_3d_xz_curv_displ_matrix(e_y, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix in the local `x-z` plane (i.e.
    bending about the `y` axis).

    The curvature $d^2w/dx^2$ is computed in the local coordinate system of the
    beam as $d^2w/dx^2 = B U$. Here $B$ is the curvature-displacement matrix and
    $U$ is the displacement vector. All three displacement components and all
    three rotation components at each joint are assumed, so the matrix $B$ has
    one row and twelve columns.
    """
    d2Ndxi2 = beam_3d_xz_shape_fun_xi2(xi)
    B = zeros((1, 12))
    B[0, 0:3] = d2Ndxi2[0] * (2 / h) ** 2 * e_z
    B[0, 3:6] = (h / 2) * d2Ndxi2[1] * (2 / h) ** 2 * e_y
    B[0, 6:9] = d2Ndxi2[2] * (2 / h) ** 2 * e_z
    B[0, 9:12] = (h / 2) * d2Ndxi2[3] * (2 / h) ** 2 * e_y
    return B


def beam_3d_xy_curv_displ_matrix(e_y, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix in the local `x-y` plane (i.e.
    bending about the `z` axis).

    The curvature $d^2v/dx^2$ is computed in the local coordinate system of the
    beam as $d^2v/dx^2 = B U$. Here $B$ is the curvature-displacement matrix and
    $U$ is the displacement vector. All three displacement components and all
    three rotation components at each joint are assumed, so the matrix $B$ has
    one row and twelve columns.
    """
    d2Ndxi2 = beam_3d_xy_shape_fun_xi2(xi)
    B = zeros((1, 12))
    B[0, 0:3] = d2Ndxi2[0] * (2 / h) ** 2 * e_y
    B[0, 3:6] = (h / 2) * d2Ndxi2[1] * (2 / h) ** 2 * e_z
    B[0, 6:9] = d2Ndxi2[2] * (2 / h) ** 2 * e_y
    B[0, 9:12] = (h / 2) * d2Ndxi2[3] * (2 / h) ** 2 * e_z
    return B


def beam_3d_bending_stiffness(e_y, e_z, h, E, Iy, Iz):
    """
    Compute 3d beam stiffness matrices for bending in the planes `x-y` and `x-z`.

    Two-point Gauss quadrature is used to compute the stiffness matrices.
    """
    xiG = [-1 / sqrt(3), 1 / sqrt(3)]
    WG = [1, 1]
    Kxy = zeros((12, 12))
    for q in range(2):
        B = beam_3d_xy_curv_displ_matrix(e_y, e_z, h, xiG[q])
        Kxy += E * Iz * outer(B.T, B) * WG[q] * (h / 2)
    Kxz = zeros((12, 12))
    for q in range(2):
        B = beam_3d_xz_curv_displ_matrix(e_y, e_z, h, xiG[q])
        Kxz += E * Iy * outer(B.T, B) * WG[q] * (h / 2)
    return Kxy, Kxz


def beam_2d_curv_displ_matrix(e_z, h, xi):
    """
    Compute 2d beam curvature-displacement matrix.

    Here the curvatures is with respect to the physical coordinate measured
    along the member (local `x`).

    The curvature $d^2w/dx^2$ is computed in the local coordinate system of the
    beam as $d^2w/dx^2 = B U$. Here $B$ is the curvature-displacement matrix and
    $U$ is the displacement vector. Two displacement components and one rotation
    component at each joint are assumed, so the matrix $B$ has one row and six
    columns.
    """
    d2Ndxi2 = beam_2d_shape_fun_xi2(xi)
    B = zeros((1, 6))
    B[0, 0:2] = d2Ndxi2[0] * (2 / h) ** 2 * e_z
    B[0, 2] = (h / 2) * d2Ndxi2[1] * (2 / h) ** 2
    B[0, 3:5] = d2Ndxi2[2] * (2 / h) ** 2 * e_z
    B[0, 5] = (h / 2) * d2Ndxi2[3] * (2 / h) ** 2
    return B


def beam_2d_3rd_deriv_displ_matrix(e_z, h, xi):
    """
    Compute beam third derivative-displacement matrix.

    Here the third derivative is with respect to the physical coordinate measured
    along the member (local `x`).

    The third derivative $d^3w/dx^3$ is computed in the local coordinate system of the
    beam as $d^3w/dx^3 = B U$. Here $B$ is the third-derivative-displacement matrix and
    $U$ is the displacement vector. Two displacement components and one rotation
    component at each joint are assumed, so the matrix $B$ has one row and six
    columns.
    """
    d2Ndxi3 = beam_2d_shape_fun_xi3(xi)
    B = zeros((1, 6))
    B[0, 0:2] = d2Ndxi3[0] * (2 / h) ** 3 * e_z
    B[0, 2] = (h / 2) * d2Ndxi3[1] * (2 / h) ** 3
    B[0, 3:5] = d2Ndxi3[2] * (2 / h) ** 3 * e_z
    B[0, 5] = (h / 2) * d2Ndxi3[3] * (2 / h) ** 3
    return B


def beam_3d_xz_3rd_deriv_displ_matrix(e_y, e_z, h):
    """
    Compute 3d beam third derivative-displacement matrix for displacements in
    the `x-z` plane.

    Here the third derivative is with respect to the physical coordinate measured
    along the member (local `x`).

    The third derivative $d^3w/dx^3$ is computed in the local coordinate system of the
    beam as $d^3w/dx^3 = B U$. Here $B$ is the third-derivative-displacement matrix and
    $U$ is the displacement vector. All three displacement components and all
    three rotation components at each joint are assumed, so the matrix $B$ has
    one row and twelve columns.
    """
    d2Ndxi3 = beam_3d_xz_shape_fun_xi3(0.0)
    B = zeros((1, 12))
    B[0, 0:3] = d2Ndxi3[0] * (2 / h) ** 3 * e_z
    B[0, 3:6] = (h / 2) * d2Ndxi3[1] * (2 / h) ** 3 * e_y
    B[0, 6:9] = d2Ndxi3[2] * (2 / h) ** 3 * e_z
    B[0, 9:12] = (h / 2) * d2Ndxi3[3] * (2 / h) ** 3 * e_y
    return B


def beam_3d_xy_3rd_deriv_displ_matrix(e_y, e_z, h):
    """
    Compute 3d beam third derivative-displacement matrix for displacements in
    the `x-y` plane.

    Here the third derivative is with respect to the physical coordinate measured
    along the member (local `x`).

    The third derivative $d^3v/dx^3$ is computed in the local coordinate system of the
    beam as $d^3v/dx^3 = B U$. Here $B$ is the third-derivative-displacement matrix and
    $U$ is the displacement vector. All three displacement components and all
    three rotation components at each joint are assumed, so the matrix $B$ has
    one row and twelve columns.
    """
    d2Ndxi3 = beam_3d_xy_shape_fun_xi3(0.0)
    B = zeros((1, 12))
    B[0, 0:3] = d2Ndxi3[0] * (2 / h) ** 3 * e_y
    B[0, 3:6] = (h / 2) * d2Ndxi3[1] * (2 / h) ** 3 * e_z
    B[0, 6:9] = d2Ndxi3[2] * (2 / h) ** 3 * e_y
    B[0, 9:12] = (h / 2) * d2Ndxi3[3] * (2 / h) ** 3 * e_z
    return B


def beam_2d_moment(member, i, j, xi):
    """
    Compute 2d beam moment based on the displacements stored at the joints.
    The moment is computed at the parametric location `xi` along the beam.

    The moment is mathematically defined as $M = -EI d^2w/dx^2$.

    The curvature is computed with the curvature-displacement matrix $B$ by the function
    `beam_2d_curv_displ_matrix`.
    """
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    sect = member["section"]
    E, I = sect["E"], sect["I"]
    ui, uj = i["displacements"], j["displacements"]
    u = concatenate([ui, uj])
    B = beam_2d_curv_displ_matrix(e_z, h, xi)
    return -E * I * dot(B, u)


def beam_3d_moment(member, i, j, axis, xi):
    """
    Compute 3d beam moment based on the displacements stored at the joints. The
    moment is computed at the parametric location `xi` along the beam. The
    moment acts about the axis specified by the string `axis` (`'y'` or `'z'`).

    The moments are mathematically defined as $M_y = -EI_y d^2w/dx^2$ for bending
    about the `y` axis, and $M_z = +EI_z d^2v/dx^2$ for bending about the `z`
    axis.

    The curvatures are computed with  curvature-displacement matrices $B$ by the
    functions `beam_3d_xz_curv_displ_matrix` and `beam_3d_xy_curv_displ_matrix` respectively.
    """
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    E, Iy, Iz = sect["E"], sect["Iy"], sect["Iz"]
    ui, uj = i["displacements"], j["displacements"]
    u = concatenate([ui, uj])
    if axis == "y":
        B = beam_3d_xz_curv_displ_matrix(e_y, e_z, h, xi)
        M = -E * Iy * dot(B, u)
    else:
        B = beam_3d_xy_curv_displ_matrix(e_y, e_z, h, xi)
        M = +E * Iz * dot(B, u)
    return M


def beam_3d_torsion_moment(member, i, j):
    """
    Compute 3d beam torsion moment based on the displacements stored at the joints.
    The moment is uniform along the beam.
    """
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    G, J = sect["G"], sect["J"]
    ui, uj = i["displacements"][3:6], j["displacements"][3:6]
    u = concatenate([ui, uj])
    B = beam_3d_torsion_displ_matrix(e_x, h)
    T = G * J * dot(B, u)
    return T


def beam_2d_axial_force(member, i, j):
    """
    Compute 2d beam axial force based on the displacements stored at the joints.
    The axial force is uniform along the beam.
    """
    sect = member["section"]
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    E, A = sect["E"], sect["A"]
    ui, uj = i["displacements"][0:2], j["displacements"][0:2]
    u = concatenate([ui, uj])
    B = truss.truss_strain_displacement(e_x, h)
    N = E * A * dot(B, u)
    return N


def beam_3d_axial_force(member, i, j):
    """
    Compute 3d beam axial force based on the displacements stored at the joints.
    The axial force is uniform along the beam.
    """
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    E, A = sect["E"], sect["A"]
    ui, uj = i["displacements"][0:3], j["displacements"][0:3]
    u = concatenate([ui, uj])
    B = beam_3d_stretch_displ_matrix(e_x, h)
    N = E * A * dot(B, u)
    return N


def beam_3d_shear_force(member, i, j, axis, xi):
    """
    Compute 3d shear force based on the displacements stored at the joints. The
    shear force in the direction of axis `axis`  (`'y'` or `'z'`) is uniform along the beam.
    """
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    E, Iy, Iz = sect["E"], sect["Iy"], sect["Iz"]
    ui, uj = i["displacements"], j["displacements"]
    u = concatenate([ui, uj])
    if axis == "z":
        B = beam_3d_xz_3rd_deriv_displ_matrix(e_y, e_z, h)
        Q = -E * Iy * dot(B, u)
    else:
        B = beam_3d_xy_3rd_deriv_displ_matrix(e_y, e_z, h)
        Q = -E * Iz * dot(B, u)
    return Q


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
    B = beam_2d_3rd_deriv_displ_matrix(e_z, h, xi)
    return -E * I * dot(B, u)


def beam_3d_stretch_displ_matrix(e_x, h):
    """
    Compute beam stretch-displacement matrix.

    The job is delegated to the truss module.
    """
    return truss.truss_strain_displacement(e_x, h)


def beam_3d_axial_stiffness(e_x, h, E, A):
    """
    Compute axial stiffness matrix.
    """
    B = beam_3d_stretch_displ_matrix(e_x, h)
    return E * A * outer(B.T, B) * h


def beam_3d_torsion_stiffness(e_x, h, G, J):
    """
    Compute torsion stiffness matrix.
    """
    B = beam_3d_torsion_displ_matrix(e_x, h)
    return G * J * outer(B.T, B) * h


def beam_3d_torsion_displ_matrix(e_x, h):
    """
    Compute torsion-displacement matrix.
    """
    B = zeros((1, 6))
    B[0, 0:3] = -e_x / h
    B[0, 3:6] = e_x / h
    return B


def assemble_stiffness(Kg, member, i, j):
    """
    Assemble beam stiffness matrix.

    Beam is here considered as a superposition of four mechanisms  -- axial bar,
    torsion bar, bending in the `x-y` plane, and bending in the `x-z` plane.
    """
    # Add stiffness in bending.
    beam_is_2d = len(i["coordinates"]) == len(j["coordinates"]) == 2
    dof = concatenate([i["dof"], j["dof"]])
    if beam_is_2d:
        e_x, e_z, h = beam_2d_member_geometry(i, j)
        sect = member["section"]
        E, I = sect["E"], sect["I"]
        k = beam_2d_bending_stiffness(e_z, h, E, I)
        Kg = assemble.assemble(Kg, dof, k)
    else:
        sect = member["section"]
        e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
        E, Iy, Iz = sect["E"], sect["Iy"], sect["Iz"]
        kxy, kxz = beam_3d_bending_stiffness(e_y, e_z, h, E, Iy, Iz)
        Kg = assemble.assemble(Kg, dof, kxy)
        Kg = assemble.assemble(Kg, dof, kxz)
    # Add stiffness in the axial direction.
    E, A = sect["E"], sect["A"]
    k = beam_3d_axial_stiffness(e_x, h, E, A)
    if beam_is_2d:
        dof = concatenate([i["dof"][0:2], j["dof"][0:2]])
    else:
        dof = concatenate([i["dof"][0:3], j["dof"][0:3]])
    Kg = assemble.assemble(Kg, dof, k)
    if not beam_is_2d:
        # Add stiffness in torsion.
        e_x, h = truss.truss_member_geometry(i, j)
        sect = member["section"]
        G, J = sect["G"], sect["J"]
        k = beam_3d_torsion_stiffness(e_x, h, G, J)
        dof = concatenate([i["dof"][3:6], j["dof"][3:6]])
        Kg = assemble.assemble(Kg, dof, k)
    return Kg


def assemble_mass(Mg, member, i, j):
    """
    Assemble beam mass matrix.
    """
    beam_is_2d = len(i["coordinates"]) == len(j["coordinates"]) == 2
    dof = concatenate([i["dof"], j["dof"]])
    sect = member["section"]
    rho, A = sect["rho"], sect["A"]
    if beam_is_2d:
        e_x, e_z, h = beam_2d_member_geometry(i, j)
        I = sect["I"]
        m = beam_2d_mass(e_x, h, rho, A, I)
    else:
        e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
        Ix, Iy, Iz = sect["Ix"], sect["Iy"], sect["Iz"]
        m = beam_3d_mass(e_x, e_y, e_z, h, rho, A, Ix, Iy, Iz)
    Mg = assemble.assemble(Mg, dof, m)
    return Mg


def beam_2d_mass(e_x, h, rho, A, I):
    """
    Compute beam mass matrix.
    """
    HLIy = rho * I * h / 2.0
    n1 = len(e_x) + 1
    m = zeros((2 * n1, 2 * n1))
    for i in range(len(e_x)):
        m[i, i] = rho * A * h / 2
        m[i + n1, i + n1] = rho * A * h / 2
    m[n1 - 1, n1 - 1] = HLIy
    m[2 * n1 - 1, 2 * n1 - 1] = HLIy
    return m


def beam_3d_mass(e_x, e_y, e_z, h, rho, A, Ix, Iy, Iz):
    """
    Compute beam mass matrix.

    The matrix is lumped and rotation inertias are included.
    """
    n1 = len(e_x)
    m = zeros((4 * n1, 4 * n1))
    HLM = A * rho * h / 2.0
    HLIx = rho * Ix * h / 2.0
    HLIy = rho * Iy * h / 2.0
    HLIz = rho * Iz * h / 2.0
    m[0, 0] = HLM
    m[1, 1] = HLM
    m[2, 2] = HLM
    m[6, 6] = HLM
    m[7, 7] = HLM
    m[8, 8] = HLM
    R = zeros((3, 3))
    R[0:3, 0] = e_x
    R[0:3, 1] = e_y
    R[0:3, 2] = e_z
    msub = zeros((3, 3))
    msub[0, 0] = +HLIx
    msub[1, 1] = +HLIy
    msub[2, 2] = +HLIz
    msub = dot(dot(R, msub), R.T)
    m[3:6, 3:6] = msub
    m[9:12, 9:12] = msub
    return m


def beam_3d_end_forces(member, i, j):
    """
    Compute the end forces of a beam element in 3d.

    Dictionary with the keys 'Ni', 'Qyi', 'Qzi', 'Ti', 'Myi', 'Mzi', 'Nj',
    'Qyj', 'Qzj', 'Tj', 'Myj', 'Mzj',  is returned.
    """
    Ni = beam_3d_axial_force(member, i, j)
    Nj = Ni
    Ti = beam_3d_torsion_moment(member, i, j)
    Tj = Ti
    Myi = beam_3d_moment(member, i, j, "y", -1.0)
    Myj = beam_3d_moment(member, i, j, "y", +1.0)
    Mzi = beam_3d_moment(member, i, j, "z", -1.0)
    Mzj = beam_3d_moment(member, i, j, "z", +1.0)
    Qyi = beam_3d_shear_force(member, i, j, "y", -1.0)
    Qyj = beam_3d_shear_force(member, i, j, "y", +1.0)
    Qzi = beam_3d_shear_force(member, i, j, "z", -1.0)
    Qzj = beam_3d_shear_force(member, i, j, "z", +1.0)
    return dict(
        Ni=Ni[0],
        Qyi=Qyi[0],
        Qzi=Qzi[0],
        Ti=Ti[0],
        Myi=Myi[0],
        Mzi=Mzi[0],
        Nj=-Nj[0],
        Qyj=-Qyj[0],
        Qzj=-Qzj[0],
        Tj=-Tj[0],
        Myj=-Myj[0],
        Mzj=-Mzj[0],
    )
