"""
Define truss mechanical quantities.
"""

from math import sqrt
from numpy import reshape, outer, concatenate, zeros, dot, array
from pystran import geometry
from pystran import assemble


def truss_stiffness(e_x, h, E, A):
    """
    Compute truss stiffness matrix.

    The axial stiffness matrix is computed as $K = EA B^T B h$. Here $B$ is the
    stretch-displacement matrix, computed by `truss_strain_displacement`.
    """
    B = truss_strain_displacement(e_x, h)
    return E * A * outer(B.T, B) * h


def truss_2d_mass(e_x, e_z, h, rho, A):
    """
    Compute 2d truss mass matrix.

    The mass matrix is consistent, which means that it is computed as discrete
    form of the kinetic energy of the element, $\\int \\rho A \\dot u \\cdot
    \\dot u dx$.

    """
    xiG = [-1 / sqrt(3), 1 / sqrt(3)]
    WG = [1, 1]
    n = len(e_x) * 2
    m = zeros((n, n))
    for q in range(2):
        N = geometry.lin_basis(xiG[q])
        extN = concatenate([N[0] * e_x, N[1] * e_x])
        m += rho * A * outer(extN, extN) * WG[q] * (h / 2)

    # for i in range(n):
    #     m[i, i] = rho * A * h / 2
    return m


def truss_strain_displacement(e_x, h):
    """
    Compute truss strain-displacement matrix.

    The axial strain is computed as $\\varepsilon = B u$, using the strain
    displacement matrix $B$ and the displacement vector $u$.

    The dimension of the strain-displacement matrix depends on the number of
    space dimensions. The vector `e_x` is the unit vector along the truss
    member, so it could have one, two, or three components.
    """
    return reshape(concatenate((-e_x / h, e_x / h)), (1, 2 * len(e_x)))


def assemble_stiffness(Kg, member, i, j):
    """
    Assemble truss stiffness matrix.

    - `Kg` is the global stiffness matrix,
    - `member` is the truss member,
    - `i`, `j` are the joints.
    """
    sect = member["section"]
    E, A = sect["E"], sect["A"]
    if E <= 0.0:
        raise ValueError("Elastic modulus must be positive")
    if A <= 0.0:
        raise ValueError("Area must be positive")
    dim = len(i["coordinates"])
    if dim == 2:
        e_x, e_z, h = geometry.member_2d_geometry(i, j)
    else:
        e_x, e_y, e_z, h = geometry.member_3d_geometry(i, j, array([]))
    k = truss_stiffness(e_x, h, E, A)
    dof = concatenate([i["dof"][0:dim], j["dof"][0:dim]])
    return assemble.assemble(Kg, dof, k)


def assemble_mass(Mg, member, i, j):
    """
    Assemble truss mass matrix.
    """
    sect = member["section"]
    rho, A = sect["rho"], sect["A"]
    if rho <= 0.0:
        raise ValueError("Mass density must be positive")
    if A <= 0.0:
        raise ValueError("Area must be positive")
    dim = len(e_x)
    if dim == 2:
        e_x, e_z, h = geometry.member_2d_geometry(i, j)
    else:
        e_x, e_y, e_z, h = geometry.member_3d_geometry(i, j, array([]))
    m = truss_mass(e_x, h, rho, A)
    dim = len(e_x)
    dof = concatenate([i["dof"][0:dim], j["dof"][0:dim]])
    return assemble.assemble(Mg, dof, m)


def truss_axial_force(member, i, j):
    """
    Compute truss axial force based on the displacements stored at the joints.

    The force is computed as $N = EA B U$, where $B$ is the strain-displacement
    matrix (computed by `truss_strain_displacement`), $U$ is the displacement
    vector (so that $\\varepsilon  = BU$ is the axial strain), and $EA$ is the
    axial stiffness.
    """
    sect = member["section"]
    E, A = sect["E"], sect["A"]
    dim = len(i["coordinates"])
    if dim == 2:
        e_x, e_z, h = geometry.member_2d_geometry(i, j)
        ui, uj = i["displacements"][0:2], j["displacements"][0:2]
    else:
        e_x, e_y, e_z, h = geometry.member_3d_geometry(i, j, array([]))
        ui, uj = i["displacements"][0:3], j["displacements"][0:3]
    u = concatenate([ui, uj])
    B = truss_strain_displacement(e_x, h)
    N = E * A * dot(B, u)
    return N
