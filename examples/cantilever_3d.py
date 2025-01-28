"""
Simple cantilever under a concentrated force at the end.
"""

from context import pystran
from pystran import model
from pystran import section
from pystran import geometry
from pystran import plots
from math import sqrt

m = model.create(3)

h = 8.0
model.add_joint(m, 1, [0.0, 0.0, 0.0])
model.add_joint(m, 2, [h, 0.0, 0.0])
a = m["joints"][1]
model.add_support(a, model.U1)
model.add_support(a, model.U2)
model.add_support(a, model.U3)
model.add_support(a, model.UR1)
model.add_support(a, model.UR2)
model.add_support(a, model.UR3)
a = m["joints"][2]
model.add_support(a, model.U1)
model.add_support(a, model.U2)
model.add_support(a, model.UR1)
model.add_support(a, model.UR3)

E = 2.0e11
G = E / (2 * (1 + 0.3))
A = 6000 / 10**6
Iy = 200e6 / 10**12
Iz = Iy / 5
Ix = Iy / 5
J = 300e3 / 10**12
xz_vector = [0, 0, 1]
s1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 1, [1, 2], s1)

d = m["joints"][2]
F = 10e3
model.add_load(d, model.U3, F)

model.number_dofs(m)

model.solve_statics(m)

print([j["displacements"] for j in m["joints"].values()])

print("Reference: ", F * h**3 / (3 * E * Iy))

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 10.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
