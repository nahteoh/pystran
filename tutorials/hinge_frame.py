"""
# Example of a support-settlement problem (Section 7.4)

This example is completely solved in the book Matrix Analysis of Structures by
Robert E. Sennett, ISBN 978-1577661436. 

Displacements and internal forces are provided in the book, and we can check our
solution against these reference values.
"""

from numpy import dot
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import beam
from pystran import plots

# US customary units, inches, pounds, seconds are assumed.

# The book gives the product of the modulus of elasticity and the moment of inertia as 2.9e6.
E = 29e6
I = 100.0
A = 10.0  # cross-sectional area does not influence the results
L = 10 * 12  # span in inches

m = model.create(2)

model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [0, L])
model.add_joint(m, 5, [0, L])
model.add_joint(m, 3, [L, L])
model.add_joint(m, 4, [L, 0.0])

# The left hand side is clamped, the other joints are simply supported.
model.add_support(m["joints"][1], model.TRANSLATION_DOFS)
model.add_support(m["joints"][4], model.TRANSLATION_DOFS)


# Define the beam members.
s1 = section.beam_2d_section("s1", E, A, I)
model.add_beam_member(m, 1, [1, 2], s1)
model.add_beam_member(m, 2, [5, 3], s1)
model.add_beam_member(m, 3, [4, 3], s1)

model.add_links(m, [2, 5], model.U1)
model.add_links(m, [2, 5], model.U2)

model.add_load(m["joints"][2], model.U1, 1000.0)

model.number_dofs(m)

model.solve_statics(m)

for jid in [2, 3]:
    j = m["joints"][jid]
    print(jid, j["displacements"])

d2 = m["joints"][2]["displacements"]
d5 = m["joints"][5]["displacements"]
if norm(d2[0:2] - d5[0:2]) > 1e-3:
    raise ValueError("Incorrect displacement")

if abs(d2[0] - 0.3985) / 0.3985 > 1e-3:
    raise ValueError("Incorrect displacement")
if abs(d2[1] - 0.00041) / 0.00041 > 1e-2:
    raise ValueError("Incorrect displacement")

# member = m["beam_members"][1]
# connectivity = member["connectivity"]
# i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
# f = beam.beam_2d_end_forces(member, i, j)
# print("Member 1 end forces: ", f)
# if abs(f["Ni"]) > 1e-3:
#     raise ValueError("Incorrect force")
# if abs(f["Qzi"] / 3.9558 - 1) > 1e-3:
#     raise ValueError("Incorrect force")
# if abs(f["Myi"] / -258.92857 - 1) > 1e-3:
#     raise ValueError("Incorrect force")

# member = m["beam_members"][2]
# connectivity = member["connectivity"]
# i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
# f = beam.beam_2d_end_forces(member, i, j)
# print("Member 2 end forces: ", f)
# if abs(f["Ni"]) > 1e-3:
#     raise ValueError("Incorrect force")
# if abs(f["Qzi"] / -1.7981 - 1) > 1e-3:
#     raise ValueError("Incorrect force")
# if abs(f["Myi"] / 215.7738 - 1) > 1e-3:
#     raise ValueError("Incorrect force")

plots.plot_setup(m, set_limits=True)
plots.plot_members(m)
plots.plot_member_numbers(m)
plots.plot_joint_numbers(m)
plots.plot_beam_orientation(m, 10.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 100.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_bending_moments(m, 0.0005)
ax.set_title("Moments")
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_shear_forces(m, 0.01)
ax.set_title("Shear forces")
plots.show(m)
