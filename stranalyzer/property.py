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
