from numpy import array, zeros
import stranalyzer.property

def create():
    m = dict()
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
    return m
        
def add_truss_member(m, identifier, connectivity, properties):      
    """
    Add truss member to the model.
    """
    if (identifier in m['truss_members']):
        raise Exception("Truss member already exists") 
    else:
        m['truss_members'][identifier] = {'connectivity' : array(connectivity), 'properties' : properties}  
    return m

def add_support(j, dir, value = 0.0):
    """
    Add support to joint.
    """
    if ('supports' not in j):
        j['supports'] = dict()
    j['supports'][dir] = value
    return j
    