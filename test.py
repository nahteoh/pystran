# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 08:09:22 2024

@author: pkonl
"""

from stranalyzer import model
from stranalyzer import property

def _test_1():
    """
    >>> import stranalyzer
    >>> from stranalyzer import model
    >>> from stranalyzer import property
    >>> p = property.truss_property('steel', 200e9, 0.001)
    >>> p['name']
    'steel'
    """
    import doctest
    doctest.testmod()


def _test_2():
    """
    >>> import stranalyzer
    >>> from stranalyzer import model
    >>> from stranalyzer import property
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    {'joints': {1: {'coordinates': array([0., 0.])}}, 'truss_members': {}}
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = property.truss_property('steel', E, A)
    >>> model.add_truss_member(m, 1, [1, 2], p1)
    {'joints': {1: {'coordinates': array([0., 0.])}}, 'truss_members': {1: {'connectivity': array([1, 2]), 'properties': {'name': 'steel', 'E': 200000000000.0, 'A': 0.001}}}}
    """
    import doctest
    doctest.testmod()

def _test_3():
    """
    >>> import stranalyzer
    >>> from stranalyzer import model
    >>> from stranalyzer import property
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    {'joints': {1: {'coordinates': array([0., 0.])}}, 'truss_members': {}}
    >>> model.add_support(m['joints'][1], 0, 0.1)
    {'coordinates': array([0., 0.]), 'supports': {0: 0.1}}
    """
    import doctest
    doctest.testmod()

def _test_4():
    """
    >>> from stranalyzer import model
    >>> from stranalyzer import property
    >>> from stranalyzer import geometry
    >>> 
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    {'joints': {1: {'coordinates': array([0., 0.])}}, 'truss_members': {}}
    >>> model.add_joint(m, 2, [1.0, 1.0])
    {'joints': {1: {'coordinates': array([0., 0.])}, 2: {'coordinates': array([1., 1.])}}, 'truss_members': {}}
    >>> d = geometry.delt(m['joints'][1]['coordinates'], m['joints'][2]['coordinates'])
    >>> geometry.delt(m['joints'][1]['coordinates'], m['joints'][2]['coordinates'])
    array([1., 1.])
    """
    import doctest
    doctest.testmod()

def _test_5():
    """
    >>> import stranalyzer
    >>> from stranalyzer import model
    >>> from stranalyzer import property
    >>> from stranalyzer import geometry
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    {'joints': {1: {'coordinates': array([0., 0.])}}, 'truss_members': {}}
    >>> model.add_joint(m, 2, [1.0, 2.0])
    {'joints': {1: {'coordinates': array([0., 0.])}, 2: {'coordinates': array([1., 2.])}}, 'truss_members': {}}
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = property.truss_property('steel', E, A)
    >>> model.add_truss_member(m, 1, [1, 2], p1)
    {'joints': {1: {'coordinates': array([0., 0.])}, 2: {'coordinates': array([1., 2.])}}, 'truss_members': {1: {'connectivity': array([1, 2]), 'properties': {'name': 'steel', 'E': 200000000000.0, 'A': 0.001}}}}
    >>> i = m['truss_members'][1]['connectivity'][0]  
    >>> j = m['truss_members'][1]['connectivity'][1]  
    >>> d = geometry.delt(m['joints'][i]['coordinates'], m['joints'][j]['coordinates'])
    >>> print(d)
    [1. 2.]
    """
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test_1()
    _test_2()
    _test_3()
    _test_4()
    _test_5()

