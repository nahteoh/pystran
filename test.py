# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 08:09:22 2024

@author: pkonl
"""

from pystran import model
from pystran import section


def _test_1():
    """
    >>> import pystran
    >>> from pystran import model
    >>> from pystran import section
    >>> p = section.truss_section('steel', 200e9, 0.001)
    >>> p['name']
    'steel'
    """


def _test_2():
    """
    >>> import pystran
    >>> from pystran import model
    >>> from pystran import section
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.1, -0.2])
    >>> print(m['joints'])
    {1: {'coordinates': array([ 0.1, -0.2])}}
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = section.truss_section('steel', E, A)
    >>> model.add_truss_member(m, 1, [1, 2], p1)
    >>> print(m['truss_members'])
    {1: {'connectivity': array([1, 2], dtype=int32), 'section': {'name': 'steel', 'E': 200000000000.0, 'A': 0.001}}}
    """


def _test_3():
    """
    >>> import pystran
    >>> from pystran import model
    >>> from pystran import section
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    >>> model.add_support(m['joints'][1], 0, 0.1)
    >>> print(m['joints'][1]['supports'])
    {0: 0.1}
    >>>
    """


def _test_4():
    """
    >>> from pystran import model
    >>> from pystran import section
    >>> from pystran import geometry
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    >>> model.add_joint(m, 2, [1.0, 1.0])
    >>> d = geometry.delt(m['joints'][1]['coordinates'], m['joints'][2]['coordinates'])
    >>> geometry.delt(m['joints'][1]['coordinates'], m['joints'][2]['coordinates'])
    array([1., 1.])
    """


def _test_5():
    """
    >>> import pystran
    >>> from pystran import model
    >>> from pystran import section
    >>> from pystran import geometry
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    >>> model.add_joint(m, 2, [1.0, 2.0])
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = section.truss_section('steel', E, A)
    >>> model.add_truss_member(m, 1, [1, 2], p1)
    >>> i = m['truss_members'][1]['connectivity'][0]
    >>> j = m['truss_members'][1]['connectivity'][1]
    >>> d = geometry.delt(m['joints'][i]['coordinates'], m['joints'][j]['coordinates'])
    >>> print(d)
    [1. 2.]
    """


def _test_6():
    """
    >>> import pystran
    >>> from pystran import model
    >>> from pystran import section
    >>> from pystran import geometry
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    >>> model.add_joint(m, 2, [1.0, 2.0])
    >>> model.add_joint(m, 3, [-1.0, 1.0])
    >>> model.add_support(m['joints'][1], 0)
    >>> model.add_support(m['joints'][1], 1)
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = section.truss_section('steel', E, A)
    >>> model.add_truss_member(m, 1, [1, 2], p1)
    >>> model.add_truss_member(m, 2, [3, 2], p1)
    >>> model.number_dofs(m)
    >>> print(m['nfreedof'])
    4
    >>> print(m['ntotaldof'])
    6
    >>> print([j['dof'] for j in m['joints'].values()])
    [array([4, 5], dtype=int32), array([0, 1], dtype=int32), array([2, 3], dtype=int32)]
    >>>
    """


def _test_7():
    """
    >>> import pystran
    >>> from pystran import model
    >>> from pystran import section
    >>> from pystran import geometry
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    >>> model.add_joint(m, 2, [1.0, 2.0])
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = section.truss_section('steel', E, A)
    >>> model.add_truss_member(m, 1, [1, 2], p1)
    >>> i = m['truss_members'][1]['connectivity'][0]
    >>> j = m['truss_members'][1]['connectivity'][1]
    >>> d = geometry.delt(m['joints'][i]['coordinates'], m['joints'][j]['coordinates'])
    >>> print(d)
    [1. 2.]
    >>> L = geometry.len(m['joints'][i]['coordinates'], m['joints'][j]['coordinates'])
    >>> print(L)
    2.23606797749979
    >>>
    """


def _test_8():
    """
    >>> import pystran
    >>> from pystran import model
    >>> from pystran import section
    >>> from pystran import geometry
    >>> from numpy import array, dot, outer
    >>> m = model.create()
    >>> model.add_joint(m, 1, [0.0, 0.0])
    >>> model.add_joint(m, 2, [0.0, -2.0])
    >>> model.add_joint(m, 3, [-1.0, 0.0])
    >>> model.add_joint(m, 4, [+1.0, 0.0])
    >>> model.add_support(m['joints'][1], 0)
    >>> model.add_support(m['joints'][1], 1)
    >>> model.add_support(m['joints'][3], 0)
    >>> model.add_support(m['joints'][3], 1)
    >>> model.add_support(m['joints'][4], 0)
    >>> model.add_support(m['joints'][4], 1)
    >>> E = 200e9
    >>> A = 0.001
    >>> p1 = section.truss_section('steel', E, A)
    >>> model.add_truss_member(m, 1, [1, 2], p1)
    >>> model.add_truss_member(m, 2, [3, 2], p1)
    >>> model.add_truss_member(m, 3, [4, 2], p1)
    >>> model.number_dofs(m)
    >>> print(m['nfreedof'])
    2
    >>> print(m['ntotaldof'])
    8
    >>> print([j['dof'] for j in m['joints'].values()])
    [array([2, 3], dtype=int32), array([0, 1], dtype=int32), array([4, 5], dtype=int32), array([6, 7], dtype=int32)]
    >>>
    >>> from pystran import truss
    >>> d = geometry.delt(m['joints'][1]['coordinates'], m['joints'][2][\'coordinates'])
    >>>
    >>> B = truss.strain_displacement(d, 2.0)
    >>> print(B.shape)
    (1, 4)
    >>> K = outer(B.T, B)
    >>>
    """


if __name__ == "__main__":
    import doctest

    doctest.testmod()
