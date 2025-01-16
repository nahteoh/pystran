"""
Define utility for assembling.
"""

from numpy import arange

def assemble(Kg, dof, k):
    for r in arange(len(dof)):
        for c in arange(len(dof)):
            gr, gc = dof[r], dof[c]
            Kg[gr, gc] += k[r, c]
    return Kg
