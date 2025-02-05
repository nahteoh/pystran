"""
Define the functions for defining and manipulating degrees of freedom.
"""

from math import sqrt
import numpy
from numpy import array, zeros, dot, mean, concatenate
from numpy.linalg import norm
import scipy
import pystran.section

U1 = 0
"""
This is a designation of the degree of freedom as translation along `X`.
"""
U2 = 1
"""
This is a designation of the degree of freedom as translation along `Z` (in 2D models)
or along `Y` (in 3D models).
"""
U3 = 2
"""
This is a designation of the degree of freedom as translation along `Z` (in 3D models).
"""
UR1 = 3
"""
This is a designation of the degree of freedom as rotation about  `X` (in 3D models).
"""
UR2 = 4
"""
This is a designation of the degree of freedom as rotation about  `Y` (in 3D models).
"""
UR3 = 5
"""
This is a designation of the degree of freedom as rotation about  `Y` (in 2D
models) or rotation about `Z`  (in 3D models).
"""
ALL_DOFS = 100
"""
This is a designation of all the degrees of freedom, translations and rotations
(`U1`, `U2`, `UR3` in 2D models, or `U1`, `U2`, `U3`, `UR1`, `UR2`, `UR3` in 3D
models). It may be used to specify the clamped condition for the joint.
"""
TRANSLATION_DOFS = 200
"""
This is a designation of the translation degrees of freedom (`U1`, `U2`,  in 2D
models, or `U1`, `U2`,  `U3` in 3D models). It may be used to specify the pinned
condition for the joint.
"""


def translation_dofs(dim):
    """
    Return the translation degrees of freedom.

    The list varies according to whether `dim` implies two dimensions (2) or
    three dimensions (3): `[U1, U2]` or `[U1, U2, U3]`.
    """
    if dim == 2:
        return [U1, U2]
    else:
        return [U1, U2, U3]


def rotation_dofs(dim):
    """
    Return the rotation degrees of freedom.

    The list varies according to whether `dim` implies two dimensions (2) or
    three dimensions (3): `[UR3]` or `[UR1, UR2, UR3]`.
    """
    if dim == 2:
        return [UR3]
    else:
        return [UR1, UR2, UR3]


def translation_and_rotation_dofs(dim):
    """
    Return the translation and rotation degrees of freedom.

    See `translation_dofs` and `rotation_dofs`.
    """
    return translation_dofs(dim) + rotation_dofs(dim)


def prescribed_dofs_and_values(dim, dof, value):
    """
    Compute prescribed degrees of freedom and values for a particular support type.

    For a single `dof` and `value`, return just the tuple of the `[dof]` and
    `[value]`.

    For `dof` equal to `ALL_DOFS`, return the translation and rotation degrees
    of freedom and zero values.

    For `dof` equal to `TRANSLATION_DOFS`, return the translation degrees of
    freedom and zero values.
    """
    if dof == ALL_DOFS:
        return translation_and_rotation_dofs(dim), [
            0.0 for d in translation_and_rotation_dofs(dim)
        ]
    elif dof == TRANSLATION_DOFS:
        return translation_dofs(dim), [0.0 for d in translation_dofs(dim)]
    return [dof], [value]
