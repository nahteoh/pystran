"""
Simple geometry utilities.
"""

from math import sqrt
from numpy import array, dot
from numpy.linalg import norm, cross


def delt(xi, xj):
    """
    Compute oriented vector from the first to the second joint
    """
    return xj - xi


def vlen(xi, xj):
    """
    Compute distance from the first to the second joint
    """
    return norm(delt(xi, xj))


def lin_basis(xi):
    """
    Compute linear basis functions for an interval $-1\\le\\xi\\le+1$.
    """
    return array([(xi - 1) / (-1 - 1), (xi - -1) / (1 - -1)])


def interpolate(xi, x1, x2):
    """
    Interpolate linearly between two quantities.
    """
    N = lin_basis(xi)
    return N[0] * x1 + N[1] * x2


def herm_basis(xi):
    """
    Compute the Hermite basis functions.

    An array of basis function values is returned (i.e. $[N_1(\\xi), ..., N_4(\\xi)]$).
    """
    return array(
        [
            (2 - 3 * xi + xi**3) / 4,
            (-1 + xi + xi**2 - xi**3) / 4,
            (2 + 3 * xi - xi**3) / 4,
            (+1 + xi - xi**2 - xi**3) / 4,
        ]
    )


def herm_basis_xi(xi):
    """
    Compute the first derivative wrt $\\xi$ of the Hermite basis functions.

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


def herm_basis_xi2(xi):
    """
    Compute the second derivative wrt $\\xi$ of the Hermite basis functions.

    An array of second derivatives of shape functions is returned (i.e.
    $[d^2N_1(\\xi)/d\\xi^2, ..., d^2N_4(\\xi)/d\\xi^2]$).
    """
    return array([(6 * xi) / 4, (2 - 6 * xi) / 4, (-6 * xi) / 4, (-2 - 6 * xi) / 4])


def herm_basis_xi3(xi):
    """
    Compute the third derivative wrt $\\xi$ of the Hermite basis functions.

    An array of third derivatives of shape functions is returned (i.e.
    $[d^3N_1(\\xi)/d\\xi^3, ..., d^3N_4(\\xi)/d\\xi^3]$).
    """
    return array([(6) / 4, (-6) / 4, (-6) / 4, (-6) / 4])


def member_2d_geometry(i, j):
    """
    Compute 2d member geometry.

    A local coordinate system is attached to the member such that the `x` axis is
    along the member axis. The deformation of the member is considered in the x-z
    plane.

    Vector `e_x` is the direction vector along the axis of the member. `e_z` is
    the direction vector perpendicular to the axis of the member. These two
    vectors form a left-handed coordinate system (consistent with the sign
    convention in the book): The deflection `w` is measured positive downwards,
    while the `x` coordinate is measured left to right. So in two dimensions
    `e_x` and `e_z` form a left-handed coordinate system. In reality, the
    complete coordinate system is right-handed, as the not-used basis vector is
    `e_y`, which points out of the plane of the screen (page).
    """
    e_x = delt(i["coordinates"], j["coordinates"])
    h = vlen(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h
    e_z = array([e_x[1], -e_x[0]])
    return e_x, e_z, h


def member_3d_geometry(i, j, xz_vector):
    """
    Compute 3d member geometry.

    A local coordinate system is attached to the member such that the `x` axis
    is along the member axis. The deformation of the member is considered in the
    `x-y` and `x-z` plane.

    Vector `e_x` is the direction vector along the axis of the member. `e_z` is
    the direction vector perpendicular to the axis of the member. These two
    vectors form a right-handed coordinate system, completed by `e_y`.

    The plane `x-z` is defined by the vector `xz_vector` and the member axis
    (i.e. `e_x`). Therefore, the vector `xz_vector` must not be parallel to the
    member axis.

    - `i` and `j` = the two joints of the member.
    - `xz_vector` = the vector that defines the `x-z` plane of the member-local
      coordinate system. It does not need to be of unit length, but it must not
      be parallel to the member axis. This vector is not defined for a truss
      member, and will be passed in as empty. Heuristics will be then used to
      orient the planes.
    """
    e_x = delt(i["coordinates"], j["coordinates"])
    h = vlen(i["coordinates"], j["coordinates"])
    if h <= 0.0:
        raise ZeroDivisionError("Length of element must be positive")
    e_x /= h  # normalize the unit length
    if len(xz_vector) == 0:
        xz_vector = array([1.0, 0.0, 0.0])
        if abs(dot(e_x, xz_vector)) > 0.99 * norm(xz_vector):
            xz_vector = array([0.0, 1.0, 0.0])
    if abs(dot(e_x, xz_vector)) > 0.99 * norm(xz_vector):
        raise ZeroDivisionError("xz_vector must not be parallel to the beam axis")
    e_y = cross(xz_vector, e_x)
    e_y = e_y / norm(e_y)
    e_z = cross(e_x, e_y)
    return e_x, e_y, e_z, h
