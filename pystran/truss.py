"""
Define truss mechanical quantities.
"""

from numpy import (
    reshape,
    outer,
    concatenate,
    zeros,
)
from pystran import geometry
from pystran import assemble


def truss_member_geometry(i, j):
    """
    Compute truss geometry.
    """
    e_x = geometry.delt(i["coordinates"], j["coordinates"])
    h = geometry.vlen(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h
    return e_x, h


def stiffness(e_x, h, E, A):
    """
    Compute truss stiffness matrix.
    """
    B = strain_displacement(e_x, h)
    return E * A * outer(B.T, B) * h


def mass(e_x, h, rho, A):
    """
    Compute truss mass matrix.
    """
    n = len(e_x) * 2
    m = zeros((n, n))
    for i in range(n):
        m[i, i] = rho * A * h / 2
    return m


def strain_displacement(e_x, h):
    """
    Compute truss strain displacement matrix.
    """
    return reshape(concatenate((-e_x / h, e_x / h)), (1, 2 * len(e_x)))


def assemble_stiffness(Kg, member, i, j):
    """
    Assemble truss stiffness matrix.
    """
    e_x, h = truss_member_geometry(i, j)
    sect = member["section"]
    E, A = sect["E"], sect["A"]
    if E <= 0.0:
        raise ValueError("Elastic modulus must be positive")
    if A <= 0.0:
        raise ValueError("Area must be positive")
    k = stiffness(e_x, h, E, A)
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
