# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 08:09:22 2024

@author: pkonl
"""

from stranalyzer import model
from stranalyzer import property
from stranalyzer import geometry

m = model.create()
model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [1.0, 1.0])
print(m)

d = geometry.delt(m['joints'][1]['coordinates'], m['joints'][2]['coordinates'])
print(d)

import stranalyzer
from stranalyzer import model
from stranalyzer import property
m = model.create()
model.add_joint(m, 1, [0.0, 0.0])
E = 200e9
A = 0.001
p1 = property.truss_property('steel', E, A)
model.add_truss_member(m, 1, [1, 2], p1)


import stranalyzer
from stranalyzer import model
from stranalyzer import property
m = model.create()
model.add_joint(m, 1, [0.0, 0.0])
model.add_support(m['joints'][1], 0, 0.1)

import stranalyzer
from stranalyzer import model
from stranalyzer import property
from stranalyzer import geometry
m = model.create()
model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [1.0, 2.0])
E = 200e9
A = 0.001
p1 = property.truss_property('steel', E, A)
model.add_truss_member(m, 1, [1, 2], p1)
i = m['truss_members'][1]['connectivity'][0]  
j = m['truss_members'][1]['connectivity'][1]  
d = geometry.delt(m['joints'][i]['coordinates'], m['joints'][j]['coordinates'])
print(d)

