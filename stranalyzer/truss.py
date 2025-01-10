"""
Define truss mechanical quantities.
"""

from stranalyzer import geometry
from numpy import array, dot, reshape, transpose, hstack, vstack, arange, outer, concatenate
from numpy.linalg import norm

def stiffness(e_x, L, E, A):
    """
    Compute truss stiffness matrix.
    """
    if abs(norm(e_x) - 1.0) > 1e-6:
        raise Exception("Direction vector must be a unit vector")
    B = strain_displacement(e_x, L)
    return  E*A*outer(B.T, B)*L
    
def strain_displacement(e_x, L):
    """
    Compute truss strain displacement matrix.
    """
    return reshape(hstack([-e_x / L, e_x / L]), (2*len(e_x), 1))

def assemble_stiffness(Kg, member, i, j):
    """
    Assemble truss stiffness matrix.
    """
    e_x = geometry.delt(i['coordinates'], j['coordinates'])
    L = geometry.len(i['coordinates'], j['coordinates'])
    e_x /= L
    properties = member['properties']
    E, A = properties['E'], properties['A']
    k = stiffness(e_x, L, E, A)
    print(k)
    dof = concatenate([i['dof'], j['dof']])
    print(dof)
    for r in arange(len(dof)):
        for c in arange(len(dof)):
            gr, gc = dof[r], dof[c]
            Kg[gr, gc] += k[r, c]
    return Kg