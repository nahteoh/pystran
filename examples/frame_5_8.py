"""
Created on 01/12/2025

Example 5.8 from Matrix Structural Analysis: Second Edition 2nd Edition by
William McGuire, Richard H. Gallagher, Ronald D. Ziemian 

The section properties are not completely defined in the book.  They are
taken from example 4.8, which does not provide both second moments of area.
They are taken here as both the same.
"""

from math import sqrt
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import plots

m = model.create(3)

model.add_joint(m, 1, [0.0, 0.0, 0.0])
model.add_joint(m, 2, [5.0, 0.0, 0.0])
model.add_joint(m, 3, [5.0, 8.0, 0.0])
a = m["joints"][1]
b = m["joints"][2]
c = m["joints"][3]
model.add_support(a, freedoms.ALL_DOFS)
model.add_support(c, freedoms.ALL_DOFS)

E = 2.0e11
G = E / (2 * (1 + 0.3))
A = 4e3 / 1e6
Iy = 50e6 / 1e12
Iz = 50e6 / 1e12
Ix = Iy + Iz
J = 100e3 / 1e12
xz_vector = [0, 0, 1]
s1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 1, [1, 2], s1)
E = 2.0e11
G = E / (2 * (1 + 0.3))
A = 6e3 / 1e6
Iy = 200e6 / 1e12
Iz = 200e6 / 1e12
Ix = Iy + Iz
J = 300e3 / 1e12
xz_vector = [0, 0, 1]
xz_vector = [0, 0, 1]
s2 = section.beam_3d_section(
    "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 2, [2, 3], s2)

model.add_load(b, freedoms.U3, -5e3 - 7.5e3)
model.add_load(b, freedoms.UR2, -6.25e3)

model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

print([j["dof"] for j in m["joints"].values()])

model.solve_statics(m)

print([j["displacements"] for j in m["joints"].values()])

# print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

print(m["U"][0 : m["nfreedof"]])


if (
    norm(b["displacements"] - [0.0, 0.0, -0.02238452, 0.00419677, 0.00593197, 0.0])
    > 1.0e-5
):
    raise ValueError("Displacement calculation error")

print("Displacement calculation OK")

# print('Reference: ', [-0.02238452,  0.00419677,  0.00593197])

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 100.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
