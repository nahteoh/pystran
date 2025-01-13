# -*- coding: utf-8 -*-
"""
Created on 01/12/2025

@author: pkonl
"""

import stranalyzer
from stranalyzer import model
from stranalyzer import property
from stranalyzer import geometry
from numpy import array, dot, outer

m = model.create(2)

model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [5.0, 0.0])
model.add_joint(m, 3, [12.0, 0.0])
model.add_support(m['joints'][1], 0)
model.add_support(m['joints'][1], 1)
model.add_support(m['joints'][2], 0)
model.add_support(m['joints'][2], 1)
model.add_support(m['joints'][3], 0)
model.add_support(m['joints'][3], 1)

E = 3e10
A = 0.001
I = 1.44e-5
p1 = property.beam_2d_property('material_1', E, A, I)
model.add_beam_member(m, 1, [1, 2], p1)
E = 2.06e11
A = 0.001
I = 1.152e-5
p2 = property.beam_2d_property('material_2', E, A, I)
model.add_beam_member(m, 2, [3, 2], p2)

model.add_load(m['joints'][1], 2, 15e3)
model.add_load(m['joints'][2], 2, 25e3)
model.add_load(m['joints'][3], 2, -35e3)

model.number_dofs(m)

print('Number of free degrees of freedom = ', m['nfreedof'])
print('Number of all degrees of freedom = ', m['ntotaldof'])

print([j['dof'] for j in m['joints'].values()])

model.solve(m)

print(m['K'][0:3, 0:3])

print(m['U'][0:3])
    
