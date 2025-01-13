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

def beam_property(name, E, G, A, I1, I2, I3, J, x1x2_vector):
    """
    Define beam property.
    
    - `E`, `G`= Young's and shear modulus,
    - `A`= cross-sectional area,
    - `I2`, `I3`= central moment of inertia of the cross-section about the x2 and x3
    coordinate axis,
    - `J`= St Venant torsion constant.
    """
    p = dict()
    p['name'] = name
    p['E'] = E
    p['G'] = G
    p['A'] = A
    p['I1'] = I1
    p['I2'] = I2
    p['I3'] = I3
    p['J'] = J
    p['x1x2_vector'] = x1x2_vector
    return p
