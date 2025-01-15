"""
Define beam mechanical quantities.
"""

from stranalyzer import geometry
from stranalyzer import truss
from numpy import array, dot, reshape, transpose, hstack, vstack, arange, outer, concatenate, zeros
from numpy.linalg import norm, cross
from math import sqrt

def beam_shape_functions_2d(xi, h):
    return array([(2 - 3*xi + xi**3)/4,
                  (-1 + xi + xi**2 - xi**3)/4, 
                  (2 + 3*xi - xi**3)/4,
                  (+1 + xi - xi**2 - xi**3)/4])
    
def beam_shape_functions_xz(xi, h):
    return array([(2 - 3*xi + xi**3)/4,
                  (-1 + xi + xi**2 - xi**3)/4, 
                  (2 + 3*xi - xi**3)/4,
                  (+1 + xi - xi**2 - xi**3)/4])
    
def beam_shape_functions_xy(xi, h):
    return array([(2 - 3*xi + xi**3)/4,
                  -(-1 + xi + xi**2 - xi**3)/4, 
                  (2 + 3*xi - xi**3)/4,
                  -(+1 + xi - xi**2 - xi**3)/4])
    
def beam_2d_member_geometry(i, j):
    """
    Compute beam geometry.
    
    The deformation of the beam is considered in the x-z plane. `e_x` is the
    direction vector along the axis of the beam. `e_z` is the direction vector
    perpendicular to the axis of the beam. These two vectors form a right-handed
    coordinate system, considering `e_y` points out of the whiteboard
    (consistent with the sign convention in the book).
    """
    e_x = geometry.delt(i['coordinates'], j['coordinates'])
    h = geometry.len(i['coordinates'], j['coordinates'])
    e_x /= h
    # The orientation here reflects the sign convention in the book.
    # The deflection is measured positive downwards, while the x coordinate is measured left to right.
    # So in two dimensions e_x and e_z form a left-handed coordinate system.
    e_z = array([e_x[1], -e_x[0]])
    return e_x, e_z, h

def beam_3d_member_geometry(i, j, xz_vector):
    """
    Compute beam geometry.
    """
    e_x = geometry.delt(i['coordinates'], j['coordinates'])
    h = geometry.len(i['coordinates'], j['coordinates'])
    e_x /= h
    # The orientation here reflects the sign convention in the book.
    # The deflection is measured positive downwards, while the x coordinate is measured left to right.
    # So in two dimensions e_x and e_z form a left-handed coordinate system.
    if abs(dot(e_x, xz_vector)) > 0.99 * norm(xz_vector):
        raise Exception("xz_vector must not be parallel to the beam axis")
    e_y = cross(xz_vector, e_x)
    e_y = e_y / norm(e_y)
    e_z = cross(e_x, e_y)
    return e_x, e_y, e_z, h

def bending_stiffness_2d(e_x, e_z, h, E, I):
    """
    Compute 2d beam stiffness matrix.
    """
    if abs(norm(e_x) - 1.0) > 1e-6:
        raise Exception("Direction vector must be a unit vector")
    xiG = [-1/sqrt(3), 1/sqrt(3)]
    WG = [1, 1]
    K = zeros((6, 6))
    for q in range(2):
        B = curvature_displacement_2d(e_x, e_z, h, xiG[q])
        K += E * I * outer(B.T, B) * WG[q] * (h/2)
    return  K
    
def curvature_displacement_xz(e_x, e_y, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix in the local x-z plane.
    """
    B = zeros((1, 12))
    B[0, 0:3] = 6*xi/h**2*e_y
    B[0, 3:6] = (1 - 3*xi)/h*e_z
    B[0, 6:9] = -6*xi/h**2*e_y
    B[0, 9:12] = -(3*xi + 1)/h*e_z
    return B

def curvature_displacement_xy(e_x, e_y, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix in the local x-y plane.
    """
    B = zeros((1, 12))
    B[0, 0:3] = 6*xi/h**2*e_z
    B[0, 3:6] = -(1 - 3*xi)/h*e_y
    B[0, 6:9] = -6*xi/h**2*e_z
    B[0, 9:12] = (3*xi + 1)/h*e_y
    return B

def bending_stiffness_3d(e_x, e_y, e_z, h, E, Iy, Iz):
    """
    Compute 2d beam stiffness matrix.
    """
    if abs(norm(e_x) - 1.0) > 1e-6:
        raise Exception("Direction vector must be a unit vector")
    xiG = [-1/sqrt(3), 1/sqrt(3)]
    WG = [1, 1]
    Kxy = zeros((12, 12))
    for q in range(2):
        B = curvature_displacement_xy(e_x, e_y, e_z, h, xiG[q])
        Kxy += E * Iz * outer(B.T, B) * WG[q] * (h/2)
    Kxz = zeros((12, 12))
    for q in range(2):
        B = curvature_displacement_xz(e_x, e_y, e_z, h, xiG[q])
        Kxz += E * Iy * outer(B.T, B) * WG[q] * (h/2)
    return  Kxy, Kxz
    
def curvature_displacement_2d(e_x, e_z, h, xi):
    """
    Compute beam curvature-displacement matrix.
    """
    B = zeros((1, 6))
    B[0, 0:2] = 6*xi/h**2*e_z
    B[0, 2] = (1 - 3*xi)/h
    B[0, 3:5] = -6*xi/h**2*e_z
    B[0, 5] = -(3*xi + 1)/h
    return B

def third_deriv_displacement_2d(e_x, e_z, h, xi):
    """
    Compute beam third derivative-displacement matrix.
    """
    B = zeros((1, 6))
    B[0, 0:2] = 6/h**2*e_z/(h/2)
    B[0, 2] = (-3)/h/(h/2)
    B[0, 3:5] = -6/h**2*e_z/(h/2)
    B[0, 5] = -(3)/h/(h/2)
    return B

def beam_2d_moment(member, i, j, xi):
    """
    Compute 2d beam moment based on the displacements stored at the joints.
    The moment is computed at the parametric location `xi` along the beam.
    """
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    properties = member['properties']
    E, I = properties['E'], properties['I']
    ui, uj = i['displacements'], j['displacements']
    u = concatenate([ui, uj])
    B = curvature_displacement_2d(e_x, e_z, h, xi)
    return E * I * dot(B, u)

def beam_2d_shear_force(member, i, j, xi):
    """
    Compute 2d beam shear force based on the displacements stored at the joints.
    The shear force is computed at the parametric location `xi` along the beam.
    """
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    properties = member['properties']
    E, I = properties['E'], properties['I']
    ui, uj = i['displacements'], j['displacements']
    u = concatenate([ui, uj])
    B = third_deriv_displacement_2d(e_x, e_z, h, xi)
    return E * I * dot(B, u)

def _bending_stiffness_2d(member, i, j):
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    properties = member['properties']
    E, I = properties['E'], properties['I']
    return bending_stiffness_2d(e_x, e_z, h, E, I)
    
def _bending_stiffness_3d(member, i, j):
    properties = member['properties']
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, properties['xz_vector'])
    E, Iy, Iz = properties['E'], properties['Iy'], properties['Iz']
    return bending_stiffness_3d(e_x, e_y, e_z, h, E, Iy, Iz)

def _stiffness_truss(member, i, j):
    e_x, L = truss.truss_member_geometry(i, j)
    properties = member['properties']
    E, A = properties['E'], properties['A']
    return truss.stiffness(e_x, L, E, A)
    
def torsion_stiffness(e_x, h, G, J):
    """
    Compute torsion stiffness matrix.
    """
    if abs(norm(e_x) - 1.0) > 1e-6:
        raise Exception("Direction vector must be a unit vector")
    B = torsion_displacement(e_x, h)
    return  G*J*outer(B.T, B)*h
    
def torsion_displacement(e_x, h):
    """
    Compute torsion-displacement matrix.
    """
    B = zeros((1, 6))
    B[0, 0:3] = -e_x / h
    B[0, 3:6] = e_x / h
    return B

def _stiffness_torsion(member, i, j):
    e_x, L = truss.truss_member_geometry(i, j)
    properties = member['properties']
    E, A = properties['E'], properties['A']
    return truss.stiffness(e_x, L, E, A)
    
def assemble_stiffness(Kg, member, i, j):
    """
    Assemble beam stiffness matrix.
    """
    # Add stiffness in bending.
    beam_is_2d = len(i['coordinates']) == len(j['coordinates']) == 2
    dof = concatenate([i['dof'], j['dof']])
    if beam_is_2d:
        k = _bending_stiffness_2d(member, i, j)
        for r in arange(len(dof)):
            for c in arange(len(dof)):
                gr, gc = dof[r], dof[c]
                Kg[gr, gc] += k[r, c]
    else:
        kxy, kxz = _bending_stiffness_3d(member, i, j)
        for r in arange(len(dof)):
            for c in arange(len(dof)):
                gr, gc = dof[r], dof[c]
                Kg[gr, gc] += kxy[r, c]
                Kg[gr, gc] += kxz[r, c]
    # Add stiffness in the axial direction.
    k = _stiffness_truss(member, i, j)
    if beam_is_2d:
        dof = concatenate([i['dof'][0:2], j['dof'][0:2]])
    else:
        dof = concatenate([i['dof'][0:3], j['dof'][0:3]])
    for r in arange(len(dof)):
        for c in arange(len(dof)):
            gr, gc = dof[r], dof[c]
            Kg[gr, gc] += k[r, c]
    if not beam_is_2d:
        # Add stiffness in torsion.
        k = _stiffness_torsion(member, i, j)
        dof = concatenate([i['dof'][3:6], j['dof'][3:6]])
        for r in arange(len(dof)):
            for c in arange(len(dof)):
                gr, gc = dof[r], dof[c]
                Kg[gr, gc] += k[r, c]
    return Kg