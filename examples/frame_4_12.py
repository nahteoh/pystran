"""
Created on 01/12/2025

Continuous beam with an eccentrically loaded short bracket.

Example 4.12 from Matrix Structural Analysis: Second Edition 2nd Edition by
William McGuire, Richard H. Gallagher, Ronald D. Ziemian 

Not all section properties are provided. Assuming that the second moments of
area are the same four y and z here.
"""

from context import pystran
from pystran import model
from pystran import section
from pystran import plots
from pystran import beam

E = 200000  # SI units with lengths in millimeters
G = E / (2 * (1 + 0.3))

m = model.create(3)

model.add_joint(m, 1, [0.0, 0.0, 0.0])
model.add_joint(m, 2, [8000.0, 0.0, 0.0])
model.add_joint(m, 3, [13000.0, 0.0, 0.0])
model.add_joint(m, 4, [8000.0, 0.0, 40])

a = m["joints"][1]
model.add_support(a, model.U1)
model.add_support(a, model.U2)
model.add_support(a, model.U3)
model.add_support(a, model.UR1)
model.add_support(a, model.UR2)
model.add_support(a, model.UR3)
c = m["joints"][3]
model.add_support(c, model.U1)
model.add_support(c, model.U2)
model.add_support(c, model.U3)
model.add_support(c, model.UR1)
model.add_support(c, model.UR2)
model.add_support(c, model.UR3)

A = 6000
Iy = 200e6
Iz = Iy
Ix = Iy + Iz
J = 300e3
xz_vector = [0, 0, 1]
s1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 1, [1, 2], s1)

A = 4000
Iy = 50e6
Iz = Iy
Ix = Iy + Iz
J = 100e3
xz_vector = [0, 0, 1]
s2 = section.beam_3d_section(
    "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 2, [3, 2], s2)

# artificially increased cross section properties for the short bracket
A = 4000
Iy = 5000e6
Iz = Iy
Ix = Iy + Iz
J = 10000e3
xz_vector = [0, 1, 0]
s3 = section.beam_3d_section(
    "sect_3", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 3, [4, 2], s3)

d = m["joints"][4]
model.add_load(d, model.U2, -1e3)

# plots.plot_setup(m)
# plots.plot_members(m)
# plots.plot_beam_orientation(m, 1.0)
# plots.plot_member_numbers(m)
# plots.show(m)

model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

print([j["dof"] for j in m["joints"].values()])

model.solve_statics(m)

print([j["displacements"] for j in m["joints"].values()])

print(abs(m["joints"][2]["displacements"][1] - -0.545) / (0.545))
if abs(m["joints"][2]["displacements"][1] - -0.545) / (0.545) > 5.0e-3:
    raise ValueError("Displacement calculation error")
if abs(m["joints"][2]["displacements"][5] - -0.263e-4) / (0.263e-4) > 5.0e-3:
    raise ValueError("Displacement calculation error")


for k in m["beam_members"].keys():
    member = m["beam_members"][k]
    connectivity = member["connectivity"]
    i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
    f = beam.beam_3d_end_forces(member, i, j)
    print(f"Member {k}: ")
    print(
        f"   Joint {connectivity[0]}: N={f['Ni']:.5}, Qy={f['Qyi']:.5}, Qz={f['Qzi']:.5}, T={f['Ti']:.5}, My={f['Myi']:.5}, Mz={f['Mzi']:.5}: "
    )
    print(
        f"   Joint {connectivity[1]}: N={f['Nj']:.5}, Qy={f['Qyj']:.5}, Qz={f['Qzj']:.5}, T={f['Tj']:.5}, My={f['Myj']:.5}, Mz={f['Mzj']:.5}: "
    )

model.statics_reactions(m)

for jid in [1, 3]:
    j = m["joints"][jid]
    print(f"Joint {jid}:")
    print(
        f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, Rz={j['reactions'][2]:.5}, Mx={j['reactions'][3]:.5}, My={j['reactions'][4]:.5}, Mz={j['reactions'][5]:.5}: "
    )
if abs(m["joints"][1]["reactions"][5] - 1.7834e06) / 1.7834e06 > 0.5e-3:
    raise ValueError("Reaction calculation error")
if abs(m["joints"][3]["reactions"][5] - -1.415e06) / -1.415e06 > 0.5e-3:
    raise ValueError("Reaction calculation error")
if abs(m["joints"][1]["reactions"][3] - -2.6087e04) / -2.6087e04 > 0.5e-3:
    raise ValueError("Reaction calculation error")
if abs(m["joints"][3]["reactions"][3] - -1.3913e04) / -1.3913e04 > 0.5e-3:
    raise ValueError("Reaction calculation error")

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 1000.0)
ax = plots.plot_torsion_moments(m, scale=0.04)
ax.set_title("Torsion moment")
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
plots.show(m)
