from numpy import array, zeros
import property

def create():
    m = dict()
    m['joints'] = dict()
    m['truss_members'] = dict()
    return m

def add_joint(m, identifier, coordinates):
    """
    Add joint to the model.

    >>> m = create()
    >>> add_joint(m, 1, [0.0, 0.0])
    >>> 1 in m['joints']
    True
    >>> m['joints'][1]['coordinates']
    (0.0, 0.0)
    >>> 3 in m['joints']
    False
    """
    if (identifier in m['joints']):
        raise Exception("Joint already exists")
    else:
        m['joints'][identifier] = {'coordinates' : tuple(coordinates)}
        
def add_truss_member(m, identifier, connectivity, properties):      
    """
    Add truss member to the model.

    >>> m = create()
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = property.truss_property('steel', E, A)
    >>> add_truss_member(m, 1, [1, 2], p1)
    >>> 1 in m['truss_members']
    True
    """
    if (identifier in m['truss_members']):
        raise Exception("Truss member already exists") 
    else:
        m['truss_members'][identifier] = {'connectivity' : tuple(connectivity), 'properties' : properties}     
    
def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()