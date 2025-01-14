# -*- coding: utf-8 -*-
"""
Created on 01/12/2025

Example 7.3 from Structural Mechanics. Analytical and Numerical Approaches for
Structural Analysis by Lingyi Lu, Junbo Jia, Zhuo Tang.
"""
from context import stranalyzer
from stranalyzer import model
from stranalyzer import property
from stranalyzer import geometry
from stranalyzer import plots

m = model.create(2)

model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [-5.65, -5.65])
model.add_joint(m, 4, [8.0, 0.0])
model.add_support(m['joints'][2], model.U1)
model.add_support(m['joints'][2], model.U2)
model.add_support(m['joints'][2], model.UR3)
model.add_support(m['joints'][4], model.U1)
model.add_support(m['joints'][4], model.UR3)

E = 2.06e11
A = 9.6e-3
I = 1.152e-5
p2 = property.beam_2d_property('material_2', E, A, I)
model.add_beam_member(m, 1, [1, 2], p2)
model.add_beam_member(m, 2, [4, 1], p2)

model.add_load(m['joints'][1], model.UR3, -40e3)
model.add_load(m['joints'][1], model.U2, -50e3)

model.number_dofs(m)

print('Number of free degrees of freedom = ', m['nfreedof'])
print('Number of all degrees of freedom = ', m['ntotaldof'])

print([j['dof'] for j in m['joints'].values()])

model.solve(m)

print([j['displacements'] for j in m['joints'].values()])

print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

print(m['F'][0:m['nfreedof']])

print(m['U'][0:m['nfreedof']])
    
print(m['joints'][1]['displacements'])
print('Reference: ', [0.000236558620403,  -0.000674850902417,  -0.027039378483428])


plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 10.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_moments(m, scale=1.0e-4)
plots.show(m)
