# -*- coding: utf-8 -*-
"""
Created on 01/12/2025

Example 4.13 from 
Matrix Structural Analysis: Second Edition 2nd Edition
by William McGuire, Richard H. Gallagher, Ronald D. Ziemian 
"""
from context import pystran
from pystran import model
from pystran import property
from pystran import geometry
from pystran import plots
from math import sqrt

m = model.create(2)

model.add_joint(m, 1, [0.0, 5.0])
model.add_joint(m, 2, [8.0, 5.0])
model.add_joint(m, 3, [8.0, 0.0])
model.add_support(m['joints'][1], model.U1)
model.add_support(m['joints'][1], model.U2)
model.add_support(m['joints'][1], model.UR3)
model.add_support(m['joints'][3], model.U1)
model.add_support(m['joints'][3], model.U2)
model.add_support(m['joints'][3], model.UR3)

E = 2.0e11
A = 6000 / 10**6
I = 200e6 / 10**12
p1 = property.beam_2d_property('material_1', E, A, I)
model.add_beam_member(m, 1, [1, 2], p1)
E = 2.0e11
A = 4000 / 10**6
I = 50e6 / 10**12
p2 = property.beam_2d_property('material_2', E, A, I)
model.add_beam_member(m, 2, [3, 2], p2)

model.add_load(m['joints'][2], model.U1, 100e3 / sqrt(2))
model.add_load(m['joints'][2], model.U2, -100e3 / sqrt(2))
model.add_load(m['joints'][2], model.UR3, 50e3)

model.number_dofs(m)

print('Number of free degrees of freedom = ', m['nfreedof'])
print('Number of all degrees of freedom = ', m['ntotaldof'])

print([j['dof'] for j in m['joints'].values()])

model.solve(m)

print([j['displacements'] for j in m['joints'].values()])

# print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

print(m['U'][0:m['nfreedof']])


plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 100.0)
ax = plots.plot_shear_forces(m, scale=0.50e-3)
ax.set_title('Shear forces')
plots.show(m)
    


    
