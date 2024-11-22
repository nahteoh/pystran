import numpy
from numpy import array, zeros
import stranalyzer.property

def create(dim = 2):
    m = dict()
    m['dim'] = dim # Dimension of the model
    m['joints'] = dict()
    m['truss_members'] = dict()
    return m

def add_joint(m, identifier, coordinates):
    """
    Add joint to the model.
    """
    if (identifier in m['joints']):
        raise Exception("Joint already exists")
    else:
        m['joints'][identifier] = {'coordinates' : array(coordinates)}
    return None
        
def add_truss_member(m, identifier, connectivity, properties):      
    """
    Add truss member to the model.
    """
    if (identifier in m['truss_members']):
        raise Exception("Truss member already exists") 
    else:
        m['truss_members'][identifier] = {'connectivity' : array(connectivity, dtype=numpy.int32), 'properties' : properties}  
    return None

def add_support(j, dir, value = 0.0):
    """
    Add support to joint.
    """
    if ('supports' not in j):
        j['supports'] = dict()
    j['supports'][dir] = value
    return None

def number_dofs(m):
    """
    Number degrees of freedom.
    """
    n = 0
    for j in m['joints'].values():
        j['dof'] = array([i for i in range(m['dim'])], dtype=numpy.int32)
        for d in range(2):
            if ('supports' not in j or d not in j['supports']):
                j['dof'][d] = n
                n += 1
    m['nfreedof'] = n
    for j in m['joints'].values():
        for d in range(2):
            if ('supports' in j and d in j['supports']):
                j['dof'][d] = n
                n += 1
    m['ntotaldof'] = n
    return None
    