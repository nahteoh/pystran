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
    if m['joints'][identifier]['coordinates'].shape != (m['dim'], ):
        raise Exception("Coordinate dimension mismatch")
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

def add_load(j, dir, value):
    """
    Add load to joint.
    """
    if ('loads' not in j):
        j['loads'] = dict()
    j['loads'][dir] = value
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
    
def solve(m):
    nt, nf = m['ntotaldof'], m['nfreedof']
    # Assemble global stiffness matrix
    K = zeros((nt, nt))
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        stranalyzer.truss.assemble_stiffness(K, member, i, j)

    print(K[0:nf, 0:nf])

    # Apply boundary conditions
    F = zeros(m['ntotaldof'])
    for joint in m['joints'].values():
        if 'loads' in joint:
            for dir, value in joint['loads'].items():
                gr = joint['dof'][dir]
                F[gr] += value

    print(F[0:nf])
    
    U = zeros(m['ntotaldof'])
    # # Solve for displacements
    U[0:nf] = numpy.linalg.solve(K[0:nf, 0:nf], F[0:nf])

    # # Assign displacements back to joints
    for joint in m['joints'].values():
        joint['displacements'] = U[joint['dof']]

    return None