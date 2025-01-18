# -*- coding: utf-8 -*-
"""
Created on 01/12/2025

Example 5.8 from Matrix Structural Analysis: Second Edition 2nd Edition by
William McGuire, Richard H. Gallagher, Ronald D. Ziemian 

The section properties are not completely defined in the book.  They are
taken from example 4.8, which does not provide both second moments of area.
They are taken here as both the same.
"""
from context import pystran
from pystran import model
from pystran import property
from pystran import geometry
from pystran import plots
from math import sqrt
from numpy.linalg import norm

L = 6.0
E = 2.0e11
G = E / (2 * (1 + 0.3))
h = 0.4
b = 0.2
Iz = b * h ** 3 / 12
Iy = h * b ** 3 / 12
Ix = Iy + Iz
J = E * Iz / G

m = model.create(3)

A, B, C, D, F, G, H, I, O = 1, 2, 3, 4, 5, 6, 7, 8, 9
model.add_joint(m, O, [0.0, 0.0, 0.0])
model.add_joint(m, A, [+L/4-L, L/2, 0.0])
model.add_joint(m, B, [+L/4, L/2, 0.0])
model.add_joint(m, C, [+L/4+L, L/2, 0.0])
model.add_joint(m, D, [-L, 0.0, 0.0])
model.add_joint(m, F, [+L, 0.0, 0.0])
model.add_joint(m, G, [-L/4-L, -L/2, 0.0])
model.add_joint(m, H, [-L/4, -L/2, 0.0])
model.add_joint(m, I, [-L/4+L, -L/2, 0.0])

model.add_support(m['joints'][A], model.CLAMPED)
model.add_support(m['joints'][C], model.CLAMPED)
model.add_support(m['joints'][D], model.CLAMPED)
model.add_support(m['joints'][F], model.CLAMPED)
model.add_support(m['joints'][G], model.CLAMPED)
model.add_support(m['joints'][I], model.CLAMPED)

xz_vector = [0, 0, 1]
p1 = property.beam_property('property_1', E, G, A, Ix, Iy, Iz, J, xz_vector)
model.add_beam_member(m, 1, [A, B], p1)
model.add_beam_member(m, 2, [B, C], p1)
model.add_beam_member(m, 3, [D, O], p1)
model.add_beam_member(m, 4, [O, F], p1)
model.add_beam_member(m, 5, [G, H], p1)
model.add_beam_member(m, 6, [H, I], p1)
model.add_beam_member(m, 7, [H, O], p1)
model.add_beam_member(m, 8, [O, B], p1)

model.add_load(m['joints'][O], model.U3, -100e3)

model.number_dofs(m)

print('Number of free degrees of freedom = ', m['nfreedof'])
print('Number of all degrees of freedom = ', m['ntotaldof'])

print([j['dof'] for j in m['joints'].values()])

model.solve(m)

print([j['displacements'] for j in m['joints'].values()])

# print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

print(m['U'][0:m['nfreedof']])


# if norm(b['displacements'] - [ 0.,  0., -0.02238452,  0.00419677,  0.00593197,0.]) > 1.e-5:
#     raise ValueError('Displacement calculation error')
# else:
#     print('Displacement calculation OK')

# print('Reference: ', [-0.02238452,  0.00419677,  0.00593197])

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 300.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
    


    
