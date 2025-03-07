"""
Define beam mechanical quantities.
"""

from numpy import dot, outer, concatenate, zeros
from pystran import geometry
from pystran.geometry import herm_basis_xi2, herm_basis_xi3, herm_basis
from pystran import gauss
from pystran import assemble
from pystran import truss


def beam_3d_xz_shape_fun(xi):
    """
    Compute the beam shape functions for deflection in the `x-z` plane.

    An array of shape function values is returned (i.e. $[N_1(\\xi), ..., N_4(\\xi)]$).
    """
    return herm_basis(xi)


def beam_3d_xz_shape_fun_xi2(xi):
    """
    Compute the second derivative of the beam shape functions for deflection in
    the `x-z` plane.

    An array of second derivatives of shape functions is returned (i.e.
    $[d^2N_1(\\xi)/d\\xi^2, ..., d^2N_4(\\xi)/d\\xi^2]$).
    """
    return herm_basis_xi2(xi)


def beam_3d_xy_shape_fun(xi):
    """
    Compute the beam shape functions for deflection in the x-y plane.

    The signs of the shape functions that go with the rotations (i.e. the second
    and fourth) need to be reversed: An array of second derivatives of shape
    functions is returned (i.e. $[N_1(\\xi), -N_2(\\xi),
    N_3(\\xi), -N_4(\\xi)]$).
    """
    N = herm_basis(xi)
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
    d2Ndxi2 = herm_basis_xi2(xi)
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
    return herm_basis_xi3(xi)


def beam_3d_xy_shape_fun_xi3(xi):
    """
    Compute the third derivative of the beam shape functions with respect to
    `xi` for deflection in the `x-y` plane.

    The signs of the shape functions that go with the rotations (i.e. the second
    and fourth) need to be reversed: An array of third derivatives of shape
    functions is returned (i.e. $[d^3N_1(\\xi)/d\\xi^3, -d^3N_2(\\xi)/d\\xi^3,
    d^3N_3(\\xi)/d\\xi^3, -d^3N_4(\\xi)/d\\xi^3]$).
    """
    d3Ndxi3 = herm_basis_xi3(xi)
    d3Ndxi3[1] *= -1.0
    d3Ndxi3[3] *= -1.0
    return d3Ndxi3


def beam_2d_bending_stiffness(e_z, h, E, I_y):
    """
    Compute 2d beam stiffness matrix.

    Two-point Gauss quadrature is used to compute the stiffness matrix.

    The formula reads

    $K = (h/2) \\int_{-1}^{+1} EI_y B^T B  d\\xi$,

    where $B$ is the curvature-displacement matrix (computed by
    `beam_2d_curv_displ_matrix`), and $h/2$ is the Jacobian. `I_y` is the second
    moment of area about the `y` axis (which is orthogonal to the plane of
    bending).
    """
    xiG, WG = gauss.rule(2)
    K = zeros((6, 6))
    for xi, W in zip(xiG, WG):
        B = beam_2d_curv_displ_matrix(e_z, h, xi)
        K += E * I_y * outer(B.T, B) * W * (h / 2)
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
    Compute 3d beam stiffness matrices for bending in the planes `x-y` and
    `x-z`.

    Two-point Gauss quadrature is used to compute the stiffness matrices.

    The formula reads

    $K = (h/2) \\int_{-1}^{+1} EI_y B_{xz}^T B_{xz}  d\\xi$,

    for bending in the `x-z` plane, and

    $K = (h/2) \\int_{-1}^{+1} EI_z B_{xy}^T B_{xy}  d\\xi$,

    where $B_{xz}$ is the curvature-displacement matrix for bending in the `x-z`
    plane (computed by `beam_3d_xz_curv_displ_matrix`), $B_{xy}$ is the
    curvature-displacement matrix for bending in the `x-y` plane (computed by
    `beam_3d_xy_curv_displ_matrix`), and $h/2$ is the Jacobian. `I_y` is the
    second moment of area about the `y` axis, and `I_z`  is the second moment of
    area about the `z` axis.
    """
    xiG, WG = gauss.rule(2)
    Kxy = zeros((12, 12))
    for xi, W in zip(xiG, WG):
        B = beam_3d_xy_curv_displ_matrix(e_y, e_z, h, xi)
        Kxy += E * Iz * outer(B.T, B) * W * (h / 2)
    Kxz = zeros((12, 12))
    for xi, W in zip(xiG, WG):
        B = beam_3d_xz_curv_displ_matrix(e_y, e_z, h, xi)
        Kxz += E * Iy * outer(B.T, B) * W * (h / 2)
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
    d2Ndxi2 = herm_basis_xi2(xi)
    B = zeros((1, 6))
    B[0, 0:2] = d2Ndxi2[0] * (2 / h) ** 2 * e_z
    B[0, 2] = (h / 2) * d2Ndxi2[1] * (2 / h) ** 2
    B[0, 3:5] = d2Ndxi2[2] * (2 / h) ** 2 * e_z
    B[0, 5] = (h / 2) * d2Ndxi2[3] * (2 / h) ** 2
    return B


def beam_2d_3rd_deriv_displ_matrix(e_z, h):
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
    d2Ndxi3 = herm_basis_xi3(0.0)
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
    _, e_z, h = geometry.member_2d_geometry(i, j)
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

    The moments are mathematically defined as $M_y = -EI_y d^2w/dx^2$ for
    bending about the `y` axis, and $M_z = +EI_z d^2v/dx^2$ for bending about
    the `z` axis.

    The curvatures are computed with  curvature-displacement matrices $B$ by the
    functions `beam_3d_xz_curv_displ_matrix` and `beam_3d_xy_curv_displ_matrix`,
    respectively.
    """
    sect = member["section"]
    _, e_y, e_z, h = geometry.member_3d_geometry(i, j, sect["xz_vector"])
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
    Compute 3d beam torsion moment based on the displacements stored at the
    joints. The moment is uniform along the beam.

    The moment is mathematically defined as $T = GJ d\\theta_x/dx$. The rate of
    change of the axial rotation, $d\\theta_x/dx$, is computed with the
    torsion-displacement matrix $B$, obtained by the function
    `beam_3d_torsion_displ_matrix`.

    The torsion moment is uniform along the beam.
    """
    sect = member["section"]
    e_x, _, _, h = geometry.member_3d_geometry(i, j, sect["xz_vector"])
    G, J = sect["G"], sect["J"]
    ui, uj = i["displacements"][3:6], j["displacements"][3:6]
    u = concatenate([ui, uj])
    B = beam_3d_torsion_displ_matrix(e_x, h)
    T = G * J * dot(B, u)
    return T


def beam_2d_axial_force(member, i, j):
    """
    Compute 2d beam axial force based on the displacements stored at the joints.

    Refer to the function `truss.truss_strain_displacement` that computes the
    strain-displacement matrix for a truss member.

    The axial force is uniform along the beam.
    """
    sect = member["section"]
    e_x, _, h = geometry.member_2d_geometry(i, j)
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
    e_x, _, _, h = geometry.member_3d_geometry(i, j, sect["xz_vector"])
    E, A = sect["E"], sect["A"]
    ui, uj = i["displacements"][0:3], j["displacements"][0:3]
    u = concatenate([ui, uj])
    B = beam_3d_stretch_displ_matrix(e_x, h)
    N = E * A * dot(B, u)
    return N


def beam_3d_shear_force(member, i, j, axis, xi):
    """
    Compute 3d shear force based on the displacements stored at the joints.

    The shear force in the direction of axis `axis`  (`'y'` or `'z'`) is
    uniform along the beam. Therefore, `xi` value is immaterial.
    """
    sect = member["section"]
    _, e_y, e_z, h = geometry.member_3d_geometry(i, j, sect["xz_vector"])
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
    Compute 2d beam shear force based on the displacements stored at the
    joints.

    The shear force is computed at the parametric location `xi` along the beam.
    Therefore, `xi` value is immaterial.
    """
    _, e_z, h = geometry.member_2d_geometry(i, j)
    sect = member["section"]
    E, I = sect["E"], sect["I"]
    ui, uj = i["displacements"], j["displacements"]
    u = concatenate([ui, uj])
    B = beam_2d_3rd_deriv_displ_matrix(e_z, h)
    return -E * I * dot(B, u)


def beam_3d_stretch_displ_matrix(e_x, h):
    """
    Compute beam stretch-displacement matrix.

    Stretch here means axial strain.

    The job is delegated to the truss module.
    """
    return truss.truss_strain_displacement(e_x, h)


def beam_2d_3d_axial_stiffness(e_x, h, E, A):
    """
    Compute axial stiffness matrix.

    This function works both for 2D and 3D beams.

    The axial stiffness matrix is computed as $K = EA B^T B h$. Here $B$ is the
    stretch-displacement matrix, computed by `beam_3d_stretch_displ_matrix`.
    """
    B = beam_3d_stretch_displ_matrix(e_x, h)
    return E * A * outer(B.T, B) * h


def beam_3d_torsion_stiffness(e_x, h, G, J):
    """
    Compute torsion stiffness matrix.

    The torsion stiffness matrix is computed as $K = GJ B^T B h$. Here $B$ is
    the torsion-displacement matrix, computed by `beam_3d_torsion_displ_matrix`.
    """
    B = beam_3d_torsion_displ_matrix(e_x, h)
    return G * J * outer(B.T, B) * h


def beam_3d_torsion_displ_matrix(e_x, h):
    """
    Compute torsion-displacement matrix.

    The torsion-displacement matrix is constant.
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
        sect = member["section"]
        e_x, e_z, h = geometry.member_2d_geometry(i, j)
        # Add stiffness in bending.
        E, I = sect["E"], sect["I"]
        k = beam_2d_bending_stiffness(e_z, h, E, I)
        Kg = assemble.assemble(Kg, dof, k)
        # Add stiffness in the axial direction.
        E, A = sect["E"], sect["A"]
        k = beam_2d_3d_axial_stiffness(e_x, h, E, A)
        dof = concatenate([i["dof"][0:2], j["dof"][0:2]])
        Kg = assemble.assemble(Kg, dof, k)
    else:
        sect = member["section"]
        e_x, e_y, e_z, h = geometry.member_3d_geometry(i, j, sect["xz_vector"])
        # Add stiffness in bending.
        E, Iy, Iz = sect["E"], sect["Iy"], sect["Iz"]
        kxy, kxz = beam_3d_bending_stiffness(e_y, e_z, h, E, Iy, Iz)
        Kg = assemble.assemble(Kg, dof, kxy)
        Kg = assemble.assemble(Kg, dof, kxz)
        # Add stiffness in the axial direction.
        E, A = sect["E"], sect["A"]
        k = beam_2d_3d_axial_stiffness(e_x, h, E, A)
        dof = concatenate([i["dof"][0:3], j["dof"][0:3]])
        Kg = assemble.assemble(Kg, dof, k)
        # Add stiffness in torsion.
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
        e_x, e_z, h = geometry.member_2d_geometry(i, j)
        m = beam_2d_mass(e_x, e_z, h, rho, A)
    else:
        e_x, e_y, e_z, h = geometry.member_3d_geometry(i, j, sect["xz_vector"])
        Ix, Iy, Iz = sect["Ix"], sect["Iy"], sect["Iz"]
        m = beam_3d_mass(e_x, e_y, e_z, h, rho, A, Ix, Iy, Iz)
    Mg = assemble.assemble(Mg, dof, m)
    return Mg


def beam_2d_mass(e_x, e_z, h, rho, A):
    """
    Compute beam mass matrix.

    The mass matrix is consistent, which means that it is computed as discrete
    form of the kinetic energy of the element,

    $\\int \\rho A \\left(\\dot u \\cdot \\dot u +  \\dot w \\cdot \\dot
    w\\right)dx$

    where $\\dot u$ and $\\dot w$ are the velocities in the $x$ and $z$
    directions. Note that the rotational velocities of the cross sections do not
    play a role.

    The velocity $\\dot u$ is assumed to very linearly along the element, and
    the velocity $\\dot w$ is assumed to vary according to the Hermite shape
    functions.

    Gauss quadrature is used to compute the mass matrix.
    """
    xiG, WG = gauss.rule(4)
    n = (len(e_x) + 1) * 2
    m = zeros((n, n))
    for xi, W in zip(xiG, WG):
        N = geometry.lin_basis(xi)
        extN = concatenate([N[0] * e_x, [0.0], N[1] * e_x, [0.0]])
        m += rho * A * outer(extN, extN) * W * (h / 2)
        N = geometry.herm_basis(xi)
        extN = concatenate([N[0] * e_z, [(h / 2) * N[1]], N[2] * e_z, [(h / 2) * N[3]]])
        m += rho * A * outer(extN, extN) * W * (h / 2)
    return m


def beam_3d_mass(e_x, e_y, e_z, h, rho, A, Ix, Iy, Iz):
    """
    Compute beam mass matrix.

    The mass matrix is consistent, which means that it is computed as discrete
    form of the kinetic energy of the element,

    $\\int \\rho A \\left(\\dot u \\cdot \\dot u + \\dot v \\cdot \\dot v + \\dot w \\cdot \\dot
    w\\right)dx$

    where $\\dot u$, $\\dot v$, and $\\dot w$ are the velocities in the $x$, $y$, and $z$
    directions.

    The velocity $\\dot u$ is assumed to very linearly along the element, and
    the velocity $\\dot v$, $\\dot w$ is assumed to vary according to the Hermite shape
    functions. Rotations about $y$ and $z$ are ignored in the kinetic energy.

    For spinning of the beam about its axis, the kinetic energy is given by the formula

    $\\int \\rho Ix \\dot \\theta_x \\cdot \\dot \\theta_x dx$.

    Gauss quadrature is used to compute the mass matrix.
    """
    xiG, WG = gauss.rule(4)
    n = len(e_x) * 4
    m = zeros((n, n))
    for xi, W in zip(xiG, WG):
        N = geometry.lin_basis(xi)
        # Translations
        extN = concatenate([N[0] * e_x, [0.0, 0.0, 0.0], N[1] * e_x, [0.0, 0.0, 0.0]])
        m += rho * A * outer(extN, extN) * W * (h / 2)
        extN = concatenate([N[0] * e_y, [0.0, 0.0, 0.0], N[1] * e_y, [0.0, 0.0, 0.0]])
        m += rho * A * outer(extN, extN) * W * (h / 2)
        extN = concatenate([N[0] * e_z, [0.0, 0.0, 0.0], N[1] * e_z, [0.0, 0.0, 0.0]])
        m += rho * A * outer(extN, extN) * W * (h / 2)
        # Torsion
        extN = concatenate([[0.0, 0.0, 0.0], N[0] * e_x, [0.0, 0.0, 0.0], N[1] * e_x])
        m += rho * Ix * outer(extN, extN) * W * (h / 2)
        # Transverse displacements and rotations about y and z
        N = beam_3d_xz_shape_fun(xi)
        extN = concatenate(
            [N[0] * e_z, (h / 2) * N[1] * e_y, N[2] * e_z, (h / 2) * N[3] * e_y]
        )
        m += rho * Iy * outer(extN, extN) * W * (h / 2)
        N = beam_3d_xy_shape_fun(xi)
        extN = concatenate(
            [N[0] * e_y, (h / 2) * N[1] * e_z, N[2] * e_y, (h / 2) * N[3] * e_z]
        )
        m += rho * Iz * outer(extN, extN) * W * (h / 2)

    return m


def beam_3d_end_forces(member, i, j):
    """
    Compute the end forces of a beam element in 3d.

    The end forces of the beam are forces acting on the joints `i` and `j` by the beam.

    Dictionary with the keys `'Ni'`, `'Qyi'`, `'Qzi'`, `'Ti'`, `'Myi'`, `'Mzi'`, `'Nj'`,
    `'Qyj'`, `'Qzj'`, `'Tj'`, `'Myj'`, `'Mzj'`,  is returned.
    """
    Ni = beam_3d_axial_force(member, i, j)
    Nj = -Ni
    Ti = beam_3d_torsion_moment(member, i, j)
    Tj = -Ti
    Myi = beam_3d_moment(member, i, j, "y", -1.0)
    Myj = -beam_3d_moment(member, i, j, "y", +1.0)
    Mzi = beam_3d_moment(member, i, j, "z", -1.0)
    Mzj = -beam_3d_moment(member, i, j, "z", +1.0)
    Qyi = beam_3d_shear_force(member, i, j, "y", -1.0)
    Qyj = -beam_3d_shear_force(member, i, j, "y", +1.0)
    Qzi = beam_3d_shear_force(member, i, j, "z", -1.0)
    Qzj = -beam_3d_shear_force(member, i, j, "z", +1.0)
    return dict(
        Ni=Ni[0],
        Qyi=Qyi[0],
        Qzi=Qzi[0],
        Ti=Ti[0],
        Myi=Myi[0],
        Mzi=Mzi[0],
        Nj=Nj[0],
        Qyj=Qyj[0],
        Qzj=Qzj[0],
        Tj=Tj[0],
        Myj=Myj[0],
        Mzj=Mzj[0],
    )


def beam_2d_end_forces(member, i, j):
    """
     Compute the end forces of a beam element in 3d.

     The end forces of the beam are forces acting on the joints `i` and `j` by
     the beam.

     Dictionary with the keys `'Ni'`, `'Qzi'`, `'Myi'`,  `'Nj'`,
    `'Qzj'`,  `'Myj'`,  is returned.
    """
    Ni = beam_2d_axial_force(member, i, j)
    Nj = -Ni
    Myi = beam_2d_moment(member, i, j, -1.0)
    Myj = -beam_2d_moment(member, i, j, +1.0)
    Qzi = beam_2d_shear_force(member, i, j, -1.0)
    Qzj = -beam_2d_shear_force(member, i, j, +1.0)
    return dict(
        Ni=Ni[0],
        Qzi=Qzi[0],
        Myi=Myi[0],
        Nj=Nj[0],
        Qzj=Qzj[0],
        Myj=Myj[0],
    )
