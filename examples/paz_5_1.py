"""
Created on 01/12/2025

Illustrative example 5.1 from Paz, M. and Leigh, W., Structural Dynamics: Theory
and Computation, 5th ed., Springer, 2001
"""

from numpy import array, dot, pi
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import rotation
from pystran import plots


def rotx(angleindegrees):
    return rotation.rotmat3(array([+1, 0, 0]) * angleindegrees / 180.0 * pi)


L = 200.0
E = 29500.0
G = E / (2 * (1 + 0.3))

m = model.create(3)

model.add_joint(m, 1, [0.0, 0.0, 0.0])
model.add_joint(m, 12, [-L / 2, 0.0, 0.0])
model.add_joint(m, 2, [-L, 0.0, 0.0])
model.add_joint(m, 3, [0.0, -L, 0.0])
model.add_joint(m, 4, [0.0, 0.0, +L])
model.add_joint(m, 5, [0.0, 0.0, -L])

model.add_support(m["joints"][2], model.CLAMPED)
model.add_support(m["joints"][3], model.CLAMPED)
model.add_support(m["joints"][4], model.CLAMPED)
model.add_support(m["joints"][5], model.CLAMPED)

# Member 1, 2: W18x130
A, Ix, Iy, Iz, J = section.i_beam(19.3, 11.2, 1.2, 0.67)
xz_vector = dot(rotx(-30), array([0, 0, 1]))
sect_1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
# Member 1, 3: W18x130
A, Ix, Iy, Iz, J = section.i_beam(19.3, 11.2, 1.2, 0.67)
xz_vector = array([0, 0, 1])
sect_2 = section.beam_3d_section(
    "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
# Member 1, 4: W14x82
A, Ix, Iy, Iz, J = section.i_beam(14.3, 10.1, 0.855, 0.51)
xz_vector = array([1, 0, 0])
sect_3 = section.beam_3d_section(
    "sect_3", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
# Member 1, 5: W14x82
A, Ix, Iy, Iz, J = section.i_beam(14.3, 10.1, 0.855, 0.51)
xz_vector = array([0, 1, 0])
sect_4 = section.beam_3d_section(
    "sect_4", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)

model.add_beam_member(m, 1, [1, 12], sect_1)
model.add_beam_member(m, 12, [12, 2], sect_1)
model.add_beam_member(m, 2, [1, 3], sect_2)
model.add_beam_member(m, 3, [1, 4], sect_3)
model.add_beam_member(m, 4, [1, 5], sect_4)


model.add_load(m["joints"][1], model.U1, 1000)
model.add_load(m["joints"][12], model.U2, 500)
model.add_load(m["joints"][1], model.U3, -2 * L / 2)
model.add_load(m["joints"][1], model.UR1, 2 * L**2 / 12)

model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

# print([j["dof"] for j in m["joints"].values()])

model.solve(m)

# print([j["displacements"] for j in m["joints"].values()])

# print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

# print(m["U"][0 : m["nfreedof"]])# out

refsol_1 = array([0.1865, 0.0316, -0.028, 0.003095, -0.009256, -0.0275])
if norm(m["joints"][1]["displacements"] - refsol_1) > 1.0e-1 * norm(refsol_1):
    raise ValueError("Displacement calculation error")

print("Displacement calculation OK")

print("Displacements: ", m["joints"][1]["displacements"])
print("Reference: ", [0.1865, 0.0316, -0.028, 3.095e-3, -9.256e-3, -0.0275])

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_beam_orientation(m, 40)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 100.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
