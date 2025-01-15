"""
Define property dictionaries.
"""

def truss_property(name, E, A):
    """
    Define truss property.
    """
    p = dict()
    p['name'] = name
    p['E'] = E
    p['A'] = A
    return p

def beam_2d_property(name, E, A, I):
    """
    Define 2d beam property.
    
    - `E`= Young's modulus,
    - `A`= cross-sectional area,
    - `I`= central moment of inertia of the cross-section about the x3
    coordinate axis (i.e. the axis perpendicular to the plane of the beam).
    """
    p = dict()
    p['name'] = name
    p['E'] = E
    p['A'] = A
    p['I'] = I
    return p

def beam_property(name, E, G, A, Ix, Iy, Iz, J, xz_vector):
    """
    Define beam property.
    
    - `E`, `G`= Young's and shear modulus,
    - `A`= cross-sectional area,
    - `Ix`= central moment of inertia of the cross-section about the local x.
    - `Iy`, `Iz`= central moment of inertia of the cross-section about the local y and local z
    coordinate axis,
    - `J`= St Venant torsion constant.
    - `xz_vector`= vector that lies in the local x and z coordinate plane.
    """
    p = dict()
    p['name'] = name
    p['E'] = E
    p['G'] = G
    p['A'] = A
    p['Ix'] = Ix
    p['Iy'] = Iy
    p['Iz'] = Iz
    p['J'] = J
    p['xz_vector'] = xz_vector
    return p
