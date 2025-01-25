"""
Define truss mechanical quantities.
"""

from numpy import reshape, outer, concatenate, zeros, dot
from pystran import geometry
from pystran import assemble


def truss_member_geometry(i, j):
    """
    Compute truss geometry.

    Vector `e_x` is the unit vector along the truss member, and `h` is the
    length of the member.
    """
    e_x = geometry.delt(i["coordinates"], j["coordinates"])
    h = geometry.vlen(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h
    return e_x, h


def truss_stiffness(e_x, h, E, A):
    """
    Compute truss stiffness matrix.

    The axial stiffness matrix is computed as $K = EA B^T B h$. Here $B$ is the
    stretch-displacement matrix, computed by `truss_strain_displacement`.
    """
    B = truss_strain_displacement(e_x, h)
    return E * A * outer(B.T, B) * h


def truss_mass(e_x, h, rho, A):
    """
    Compute truss mass matrix.

    The mass each joint gets is computed as $m = \\rho A h / 2$.
    """
    n = len(e_x) * 2
    m = zeros((n, n))
    for i in range(n):
        m[i, i] = rho * A * h / 2
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
    e_x, h = truss_member_geometry(i, j)
    sect = member["section"]
    E, A = sect["E"], sect["A"]
    if E <= 0.0:
        raise ValueError("Elastic modulus must be positive")
    if A <= 0.0:
        raise ValueError("Area must be positive")
    k = truss_stiffness(e_x, h, E, A)
    dim = len(e_x)
    dof = concatenate([i["dof"][0:dim], j["dof"][0:dim]])
    return assemble.assemble(Kg, dof, k)


def assemble_mass(Mg, member, i, j):
    """
    Assemble truss mass matrix.
    """
    e_x, h = truss_member_geometry(i, j)
    sect = member["section"]
    rho, A = sect["rho"], sect["A"]
    if rho <= 0.0:
        raise ValueError("Mass density must be positive")
    if A <= 0.0:
        raise ValueError("Area must be positive")
    m = mass(e_x, h, rho, A)
    dim = len(e_x)
    dof = concatenate([i["dof"][0:dim], j["dof"][0:dim]])
    return assemble.assemble(Mg, dof, m)


def truss_axial_force(member, i, j):
    """
    Compute truss axial force based on the displacements stored at the joints.
    """
    sect = member["section"]
    e_x, h = truss_member_geometry(i, j)
    E, A = sect["E"], sect["A"]
    ui, uj = i["displacements"][0:3], j["displacements"][0:3]
    u = concatenate([ui, uj])
    B = truss_strain_displacement(e_x, h)
    N = E * A * dot(B, u)
    return N
