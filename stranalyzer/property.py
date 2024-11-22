"""
Define property dictionaries.
"""

def truss_property(name, E, A):
    """
    Define truss property.

    >>> p = truss_property('steel', 200e9, 0.001)
    >>> p['name']
    'steel'
    """
    p = dict()
    p['name'] = name
    p['E'] = E
    p['A'] = A
    return p


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()