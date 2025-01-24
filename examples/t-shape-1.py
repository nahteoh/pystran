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
from pystran import beam

# SI units
L = 3.0
E = 200e9
G = 80e9
A = 0.01
Iz = 1e-3
Iy = 1e-3
Ix = 2e-3
J = Ix  # Torsional constant
P = 60000

m = model.create(3)

model.add_joint(m, 1, [0.0, 0.0, 0.0])
model.add_joint(m, 2, [0.0, L, 0.0])
model.add_joint(m, 3, [0.0, L, L])
model.add_joint(m, 4, [0.0, L, -L])

model.add_support(m["joints"][1], model.CLAMPED)

xz_vector = [1, 0, 0]
s1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 1, [1, 2], s1)
model.add_beam_member(m, 2, [2, 3], s1)
model.add_beam_member(m, 3, [2, 4], s1)

model.add_load(m["joints"][3], model.U1, -P)

model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

# print([j['dof'] for j in m['joints'].values()])

model.solve_statics(m)

for id in [2, 3, 4]:
    j = m["joints"][id]
    print(id, j["displacements"])


for k in m["beam_members"].keys():
    member = m["beam_members"][k]
    connectivity = member["connectivity"]
    i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
    f = beam.beam_3d_end_forces(member, i, j)
    print(f"Member {k}: ")
    print(
        f" Joint {connectivity[0]}: N={f['Ni']:.5}, Qy={f['Qyi']:.5}, Qz={f['Qzi']:.5}, T={f['Ti']:.5}, My={f['Myi']:.5}, Mz={f['Mzi']:.5}: "
    )
    print(
        f" Joint {connectivity[1]}: N={f['Nj']:.5}, Qy={f['Qyj']:.5}, Qz={f['Qzj']:.5}, T={f['Tj']:.5}, My={f['Myj']:.5}, Mz={f['Mzj']:.5}: "
    )

model.statics_reactions(m)

for jid in [1]:
    j = m["joints"][jid]
    print(f"Joint {jid}:")
    print(
        f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, Rz={j['reactions'][2]:.5}, Mx={j['reactions'][3]:.5}, My={j['reactions'][4]:.5}, Mz={j['reactions'][5]:.5}: "
    )

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_member_numbers(m)
plots.plot_deformations(m, 100.0)
plots.plot_beam_orientation(m, 0.5)
# plots.plot_moments(m, 0.00001, "y")
# plots.plot_moments(m, 0.00001, "z")
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
