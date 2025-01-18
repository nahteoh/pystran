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
from pystran import section
from pystran import geometry
from pystran import plots
from math import sqrt
from numpy.linalg import norm

L = 96.0  # inches
E = 1e7  # psi
G = 4e6  # psi
A = 9.0  # in^2
Iz = 80.0  # in^4
Iy = 28.0  # in^4
Ix = 64.0  # in^4
J = Ix
P = 5000.0  # lb

m = model.create(3)

jA, jB, jC, jD, jE = 1, 2, 3, 4, 5
model.add_joint(m, jA, [0.0, 0.0, 0.0])
model.add_joint(m, jB, [0.0, L, 0.0])
model.add_joint(m, jC, [2 * L, L, 0.0])
model.add_joint(m, jD, [3 * L, 0.0, L])
model.add_joint(m, jE, [L, L, 0.0])

model.add_support(m["joints"][jA], model.CLAMPED)
model.add_support(m["joints"][jD], model.CLAMPED)

xz_vector = [0, 0, 1]
p1 = section.beam_3d_section("property_1", E, G, A, Ix, Iy, Iz, J, xz_vector)
model.add_beam_member(m, 1, [jA, jB], p1)
model.add_beam_member(m, 2, [jE, jB], p1)
model.add_beam_member(m, 3, [jE, jC], p1)
model.add_beam_member(m, 4, [jC, jD], p1)

model.add_load(m["joints"][jB], model.U1, 2 * P)
model.add_load(m["joints"][jE], model.U3, 4 * P)
model.add_load(m["joints"][jC], model.U2, -P)
model.add_load(m["joints"][jC], model.UR3, -P * L)

model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

# print([j['dof'] for j in m['joints'].values()])

model.solve(m)

print([j["displacements"] for j in m["joints"].values()])

# print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

# print(m['U'][0:m['nfreedof']])


# if norm(b['displacements'] - [ 0.,  0., -0.02238452,  0.00419677,  0.00593197,0.]) > 1.e-5:
#     raise ValueError('Displacement calculation error')
# else:
#     print('Displacement calculation OK')

# print('Reference: ', [-0.02238452,  0.00419677,  0.00593197])

plots.plot_setup(m)
plots.plot_members(m)
# plots.plot_member_numbers(m)
plots.plot_deformations(m, 3.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
