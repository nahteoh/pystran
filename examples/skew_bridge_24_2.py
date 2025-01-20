"""
Created on 01/12/2025

Structural Analysis: A Unified Classical and Matrix, Ghali, Amin; Neville,
Adam -- Edition 7, 2017, Taylor and Francis

Example 24.2 - Skew Bridge
"""

from context import pystran
from pystran import model
from pystran import section
from pystran import plots

L = 6.0
E = 2.0e10
G = E / (2 * (1 + 0.3))
h = 0.2
b = 0.2
Iz = b * h**3 / 12
Iy = h * b**3 / 12
Ix = Iy + Iz
J = 0.6 * E * Iy / G
q = 1.0e3

m = model.create(3)

A, B, C, D, F, G, H, I, O = 1, 2, 3, 4, 5, 6, 7, 8, 9
model.add_joint(m, O, [0.0, 0.0, 0.0])
model.add_joint(m, A, [+L / 4 - L, L / 2, 0.0])
model.add_joint(m, B, [+L / 4, L / 2, 0.0])
model.add_joint(m, C, [+L / 4 + L, L / 2, 0.0])
model.add_joint(m, D, [-L, 0.0, 0.0])
model.add_joint(m, F, [+L, 0.0, 0.0])
model.add_joint(m, G, [-L / 4 - L, -L / 2, 0.0])
model.add_joint(m, H, [-L / 4, -L / 2, 0.0])
model.add_joint(m, I, [-L / 4 + L, -L / 2, 0.0])

model.add_support(m["joints"][A], model.CLAMPED)
model.add_support(m["joints"][C], model.CLAMPED)
model.add_support(m["joints"][D], model.CLAMPED)
model.add_support(m["joints"][F], model.CLAMPED)
model.add_support(m["joints"][G], model.CLAMPED)
model.add_support(m["joints"][I], model.CLAMPED)

xz_vector = [0, 0, 1]
s1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 1, [A, B], s1)
model.add_beam_member(m, 2, [B, C], s1)
model.add_beam_member(m, 3, [D, O], s1)
model.add_beam_member(m, 4, [O, F], s1)
model.add_beam_member(m, 5, [G, H], s1)
model.add_beam_member(m, 6, [H, I], s1)
model.add_beam_member(m, 7, [H, O], s1)
model.add_beam_member(m, 8, [O, B], s1)

model.add_load(m["joints"][O], model.U3, -L * q)

model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

model.solve(m)

print(m["joints"][O]["displacements"])
print(m["joints"][B]["displacements"])

# if norm(b['displacements'] - [ 0.,  0., -0.02238452,  0.00419677,  0.00593197,0.]) > 1.e-5:
#     raise ValueError('Displacement calculation error')
# else:
#     print('Displacement calculation OK')

print("Deflection: ", 10.68 * L * (q * L**3) / (1000 * E * Iy))
print("Deflection: ", 20.31 * L * (q * L**3) / (1000 * E * Iy))
# print("Reference: ", -7.45 * (q * L**3) / (1000 * E * Iy))
# print("Reference: ", 21.13 * (q * L**3) / (1000 * E * Iy))

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 100.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
