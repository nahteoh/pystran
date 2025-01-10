# -*- coding: utf-8 -*-
"""
Analysis of Geometrically 
Nonlinear Structures 
Second Edition 

by 
Robert Levy 
Technion-Israel Institute of Technology, 
Haifa, Israel 

and 
William R. Spillers 
New Jersey Institute of Technology, 

Three bar example on page 32
"""
import stranalyzer
from stranalyzer import model
from stranalyzer import property
from stranalyzer import geometry
from numpy import array, dot, outer

m = model.create()

model.add_joint(m, 1, [10.0, 20.0])
model.add_joint(m, 2, [0.0, 20.0])
model.add_joint(m, 3, [0.0, 10.0])
model.add_joint(m, 4, [+10.0, 0.0])

E = 30000000.0
A = 0.65700000
p1 = property.truss_property('steel', E, A)
model.add_truss_member(m, 1, [1, 2], p1)
model.add_truss_member(m, 2, [1, 3], p1)
model.add_truss_member(m, 3, [1, 4], p1)

for i in [2, 3, 4]:
    for d in range(2):
        model.add_support(m['joints'][i], d)
        
model.add_load(m['joints'][1], 0, -10000.0)
model.add_load(m['joints'][1], 1, -10000.0 / 2.0)

model.number_dofs(m)
print('Total Degrees of Freedom = ', m['ntotaldof'])
print('Free Degrees of Freedom = ', m['nfreedof'])

model.solve(m)

for j in m['joints'].values():
    print(j['displacements'])
    
print('Correct solution: ', (-0.0033325938, -0.001591621))
    
