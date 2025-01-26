"""
Created on 01/12/2025

Right angle portal frame.

Example 4.13 from 
Matrix Structural Analysis: Second Edition 2nd Edition
by William McGuire, Richard H. Gallagher, Ronald D. Ziemian 
"""

from math import sqrt
from numpy import array
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import plots

m = model.create(2)

# We are working in a SI units (m)

model.add_joint(m, 1, [0.0, 5.0])
model.add_joint(m, 2, [8.0, 5.0])
model.add_joint(m, 3, [8.0, 0.0])

model.add_support(m["joints"][1], model.U1)
model.add_support(m["joints"][1], model.U2)
model.add_support(m["joints"][1], model.UR3)
model.add_support(m["joints"][3], model.U1)
model.add_support(m["joints"][3], model.U2)
model.add_support(m["joints"][3], model.UR3)

E = 2.0e11
A = 6000 / 10**6
I = 200e6 / 10**12
s1 = section.beam_2d_section("material_1", E, A, I)
model.add_beam_member(m, 1, [1, 2], s1)
E = 2.0e11
A = 4000 / 10**6
I = 50e6 / 10**12
s2 = section.beam_2d_section("material_2", E, A, I)
model.add_beam_member(m, 2, [3, 2], s2)

model.add_load(m["joints"][2], model.U1, 100e3 / sqrt(2))
model.add_load(m["joints"][2], model.U2, -100e3 / sqrt(2))
model.add_load(m["joints"][2], model.UR3, 50e3)

model.number_dofs(m)

model.solve_statics(m)

print([j["displacements"] for j in m["joints"].values()])

if (
    norm(m["joints"][2]["displacements"] - [0.00044147, -0.00039989, 0.00169432])
    / norm([0.00044147, -0.00039989, 0.00169432])
    > 1e-5
):
    raise ValueError("Displacement error")

model.statics_reactions(m)

for jid in [1, 3]:
    j = m["joints"][jid]
    print(f"Joint {jid}:")
    print(
        f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, My={j['reactions'][2]:.5} "
    )
r = m["joints"][1]["reactions"]
if (
    norm(array([r[0], r[1], r[2]]) - array([-6.6221e04, 6728.6, 1.8443e04]))
    / norm(array([-6.6221e04, 6728.6, 1.8443e04]))
    > 1e-3
):
    raise ValueError("Reaction error")
r = m["joints"][3]["reactions"]
if (
    norm(array([r[0], r[1], r[2]]) - array([-4490.2, 6.3982e04, 7836.8]))
    / norm(array([-4490.2, 6.3982e04, 7836.8]))
    > 1e-3
):
    raise ValueError("Reaction error")


plots.plot_setup(m)
plots.plot_members(m)
plots.plot_loads(m, scale=0.00001, radius=1.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 100.0)
ax = plots.plot_shear_forces(m, scale=0.50e-3)
ax.set_title("Shear forces")
plots.show(m)


plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 100.0)
ax = plots.plot_axial_forces(m, scale=0.50e-5)
ax.set_title("Axial forces")
plots.show(m)
