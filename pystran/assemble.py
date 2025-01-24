"""
Define utility for assembling.
"""

from numpy import arange


def assemble(Kg, dof, k):
    """
    Assemble local stiffness matrix `k` into global stiffness matrix `Kg`, using
    the array of degrees of freedom, `dof`, for both the rows and columns.
    In other words, `k` must be symmetric.
    """
    for r in arange(len(dof)):
        for c in arange(len(dof)):
            gr, gc = dof[r], dof[c]
            Kg[gr, gc] += k[r, c]
    return Kg
