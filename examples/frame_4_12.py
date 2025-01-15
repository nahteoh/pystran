# -*- coding: utf-8 -*-
"""
Created on 01/12/2025

Example 4.13 from 
Matrix Structural Analysis: Second Edition 2nd Edition
by William McGuire, Richard H. Gallagher, Ronald D. Ziemian 
"""
from context import stranalyzer
from stranalyzer import model
from stranalyzer import property
from stranalyzer import geometry
from stranalyzer import plots
from math import sqrt

m = model.create(3)

model.add_joint(m, 1, [0.0, 0.0, 0.0])
model.add_joint(m, 2, [8.0, 0.0, 0.0])
model.add_joint(m, 3, [13.0, 0.0, 0.0])
model.add_joint(m, 4, [8.0, 0.0, 0.04])
a = m['joints'][1]
model.add_support(a, model.U1)
model.add_support(a, model.U2)
model.add_support(a, model.U3)
model.add_support(a, model.UR1)
model.add_support(a, model.UR2)
model.add_support(a, model.UR3)
c = m['joints'][3]
model.add_support(c, model.U1)
model.add_support(c, model.U2)
model.add_support(c, model.U3)
model.add_support(c, model.UR1)
model.add_support(c, model.UR2)
model.add_support(c, model.UR3)

E = 2.0e11
G = E / (2 * (1 + 0.3))
A = 6000 / 10**6
I2 = 200e6 / 10**12
I3 = I2 / 2
I1 = I2 / 2
J = 300e3 / 10**12
xz_vector = [0, 0, 1]
p1 = property.beam_property('property_1', E, G, A, I1, I2, I3, J, xz_vector)
model.add_beam_member(m, 1, [1, 2], p1)
E = 2.0e11
A = 4000 / 10**6
I2 = 50e6 / 10**12
I3 = I2 / 2
I1 = I2 / 2
J = 100e3 / 10**12
xz_vector = [0, 0, 1]
p2 = property.beam_property('property_2', E, G, A, I1, I2, I3, J, xz_vector)
model.add_beam_member(m, 2, [3, 2], p2)
E = 2.0e11
A = 4000 / 10**6
I2 = 5000e6 / 10**12
I3 = I2 / 2
I1 = I2 / 2
J = 10000e3 / 10**12
xz_vector = [0, 1, 0]
p3 = property.beam_property('property_3', E, G, A, I1, I2, I3, J, xz_vector)
model.add_beam_member(m, 3, [4, 2], p3)

d = m['joints'][4]
model.add_load(d, model.U2, -10e3)

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
    


    
