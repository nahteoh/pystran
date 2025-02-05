"""
Created on 01/22/2025

Weaver Jr., W., Computer Programs for Structural Analysis, page 146, problem 8.
via STAAD.Pro 2023.00.03 User Manual.

These reactions are consistent with STAAD.Pro 2023.00.03 User Manual, 
assuming Poisson ratio 0.25:

Joint 3:
   Rx=-1.1041, Ry=-0.43222, Rz=0.21731, Mx=48.785, My=-17.973, Mz=96.122:        
Joint 4:
   Rx=-0.89588, Ry=1.4322, Rz=-0.21731, Mx=123.08, My=47.246, Mz=-11.72:
"""

import numpy
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import plots
from pystran import beam
from pystran import rotation

# US customary units, inches, pounds, seconds
L = 120.0
E = 30000
nu = 0.25  # Poisson's ratio assumed in STAAD documentation
G = E / (2 * (1 + nu))
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

model.add_support(m["joints"][3], freedoms.ALL_DOFS)
model.add_support(m["joints"][4], freedoms.ALL_DOFS)

xz_vector = [0, 0, 1]
sect_1 = section.beam_3d_section(
    "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
xz_vector = [0, 0, 1]
sect_2 = section.beam_3d_section(
    "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
xz_vector = rotation.rotate(m["joints"][2], m["joints"][4], [0, 1, 0], 90)
sect_3 = section.beam_3d_section(
    "sect_3", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)

model.add_beam_member(m, 1, [1, 2], sect_1)
model.add_beam_member(m, 2, [3, 1], sect_2)
model.add_beam_member(m, 3, [2, 4], sect_3)

model.add_load(m["joints"][1], freedoms.U1, F)
model.add_load(m["joints"][2], freedoms.U2, -P)
model.add_load(m["joints"][2], freedoms.UR3, -M)

model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

print([j["dof"] for j in m["joints"].values()])

model.solve_statics(m)

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

for jid in [3, 4]:
    j = m["joints"][jid]
    print(f"Joint {jid}:")
    print(
        f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, Rz={j['reactions'][2]:.5}, Mx={j['reactions'][3]:.5}, My={j['reactions'][4]:.5}, Mz={j['reactions'][5]:.5}: "
    )

allforces = model.free_body_check(m)
print("Sum of forces and moments: ", allforces)
if norm(allforces) > 1.0e-10:
    raise ValueError("Sum of forces and moments not zero")

# m["joints"][3]["reactions"] = {
#     0: -1.10,
#     1: -0.43,
#     2: 0.22,
#     3: 48.78,
#     4: -17.97,
#     5: 96.12,
# }
# m["joints"][4]["reactions"] = {
#     0: -0.90,
#     1: 1.43,
#     2: -0.22,
#     3: 123.08,
#     4: 47.25,
#     5: -11.72,
# }
# allforces = model.free_body_check(m)
# print("Sum of forces and moments: ", allforces)
# if norm(allforces) > 1.0e-10:
#     raise ValueError("Sum of forces and moments not zero")


plots.plot_setup(m)
plots.plot_members(m)
plots.plot_beam_orientation(m, 20)
ax = plots.plot_applied_moments(m, 1.0)
ax.view_init(elev=58, azim=118)
plots.show(m)
print(ax.elev, ax.azim)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_beam_orientation(m, 20)
plots.plot_applied_forces(m, 100.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_beam_orientation(m, 20)
plots.plot_deformations(m, 80.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_beam_orientation(m, 20)
ax = plots.plot_shear_forces(m, scale=50)
ax.set_title("Shear forces")
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_beam_orientation(m, 20)
ax = plots.plot_torsion_moments(m, scale=2)
ax.set_title("Torsional moments")
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_beam_orientation(m, 20)
ax = plots.plot_axial_forces(m, scale=10)
ax.set_title("Axial forces")
plots.show(m)
