"""
Created on 01/12/2025

Example 5.8 from Matrix Structural Analysis: Second Edition 2nd Edition by
William McGuire, Richard H. Gallagher, Ronald D. Ziemian 

The section properties are not completely defined in the book.  They are
taken from example 4.8, which does not provide both second moments of area.
They are taken here as both the same.
"""

from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import plots

# US customary units, inches, pounds, seconds
L = 120.0
E = 30000
G = E / (2 * (1 + 0.3))
A = 11
Iz = 56
Iy = 56
Ix = 83
J = Ix  # Torsional constant
F = 2
P = 1
M = 120

m = model.create(3)

model.add_joint(m, 3, [0.0, 0.0, 0.0])
model.add_joint(m, 1, [0.0, L, 0.0])
model.add_joint(m, 2, [2 * L, L, 0.0])
model.add_joint(m, 4, [3 * L, 0.0, L])

model.add_support(m["joints"][3], model.CLAMPED)
model.add_support(m["joints"][4], model.CLAMPED)

xz_vector = [1, 0, 0]
sect_1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 1, [3, 1], sect_1)
xz_vector = [0, 1, 0]
sect_2 = section.beam_3d_section(
    "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 2, [1, 2], sect_2)
model.add_beam_member(m, 3, [2, 4], sect_2)

model.add_load(m["joints"][1], model.U1, F)
model.add_load(m["joints"][2], model.U2, -P)
model.add_load(m["joints"][2], model.UR3, -M)

model.number_dofs(m)

# print("Number of free degrees of freedom = ", m["nfreedof"])
# print("Number of all degrees of freedom = ", m["ntotaldof"])

# print([j['dof'] for j in m['joints'].values()])

model.solve(m)

for jid in [1, 2]:
    j = m["joints"][jid]
    print(jid, j["displacements"])

# print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

# print(m['U'][0:m['nfreedof']])

# 0.22267   0.00016  -0.17182  -0.00255   0.00217  -0.00213
# 0.22202  -0.48119  -0.70161  -0.00802   0.00101  -0.00435
ref1 = [0.22267, 0.00016, -0.17182, -0.00255, 0.00217, -0.00213]
if norm(m["joints"][1]["displacements"] - ref1) > 1.0e-1 * norm(ref1):
    raise ValueError("Displacement calculation error")
else:
    print("Displacement calculation OK")
ref2 = [0.22202, -0.48119, -0.70161, -0.00802, 0.00101, -0.00435]
if norm(m["joints"][2]["displacements"] - ref2) > 1.0e-1 * norm(ref2):
    raise ValueError("Displacement calculation error")
else:
    print("Displacement calculation OK")


plots.plot_setup(m)
plots.plot_members(m)
# plots.plot_member_numbers(m)
plots.plot_deformations(m, 30.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
