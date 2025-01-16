"""
Define truss mechanical quantities.
"""

from pystran import geometry
from pystran import assemble
from numpy import array, dot, reshape, transpose, hstack, vstack, arange, outer, concatenate
from numpy.linalg import norm

def truss_member_geometry(i, j):
    """
    Compute truss geometry.
    """
    e_x = geometry.delt(i['coordinates'], j['coordinates'])
    h = geometry.len(i['coordinates'], j['coordinates'])
    if h <= 0.0:
        raise Exception("Length of element must be positive")
    e_x /= h
    return e_x, h

def stiffness(e_x, h, E, A):
    """
    Compute truss stiffness matrix.
    """
    B = strain_displacement(e_x, h)
    return  E*A*outer(B.T, B)*h
    
def strain_displacement(e_x, h):
    """
    Compute truss strain displacement matrix.
    """
    return reshape(concatenate((-e_x / h, e_x / h)), (1, 2*len(e_x)))

def assemble_stiffness(Kg, member, i, j):
    """
    Assemble truss stiffness matrix.
    """
    e_x, h = truss_member_geometry(i, j)
    properties = member['properties']
    E, A = properties['E'], properties['A']
    k = stiffness(e_x, h, E, A)
    dof = concatenate([i['dof'], j['dof']])
    return assemble.assemble(Kg, dof, k)