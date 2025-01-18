"""
Created on 01/12/2025

Example 4.13 from 
Matrix Structural Analysis: Second Edition 2nd Edition
by William McGuire, Richard H. Gallagher, Ronald D. Ziemian 
"""

from numpy import zeros, dot
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section

m = model.create(3)

# General orientation
model.add_joint(m, 1, [-1.0, 2.0, 3.0])
# model.add_joint(m, 2, [10.0, 7.0, 8.0])
model.add_joint(m, 2, [-10.0, 7.0, -8.0])
# Default orientation
# model.add_joint(m, 1, [0.0, 0.0, 0.0])
# model.add_joint(m, 2, [10.0, 0.0, 0.0])
h = norm(m["joints"][1]["coordinates"] - m["joints"][2]["coordinates"])

E = 2.0e6
G = E / (2 * (1 + 0.3))
H = 0.13
B = 0.5
A = H * B
Iy = H * B**3 / 12
Iz = H**3 * B / 12
Ix = Iy + Iz
J = Ix
xz_vector = [0, 0, 1]
s1 = section.beam_3d_section("property_1", E, G, A, Ix, Iy, Iz, J, xz_vector)
model.add_beam_member(m, 1, [1, 2], s1)

model.number_dofs(m)

nt, nf = m["ntotaldof"], m["nfreedof"]
# Assemble global stiffness matrix
K = zeros((nt, nt))
for member in m["beam_members"].values():
    connectivity = member["connectivity"]
    i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
    pystran.beam.assemble_stiffness(K, member, i, j)

K1 = K.copy()

K = zeros((nt, nt))
# Axial force
K[0, 0] = E * A / h
K[6, 6] = E * A / h
K[0, 6] = -E * A / h
K[6, 0] = -E * A / h
# Torsion
K[3, 3] = G * J / h
K[9, 9] = G * J / h
K[3, 9] = -G * J / h
K[9, 3] = -G * J / h
# Bending in xy plane
K[1, 1] = 12 * E * Iz / h**3
K[7, 7] = 12 * E * Iz / h**3
K[1, 7] = -12 * E * Iz / h**3
K[7, 1] = -12 * E * Iz / h**3
K[1, 5] = 6 * E * Iz / h**2
K[5, 1] = 6 * E * Iz / h**2
K[1, 11] = 6 * E * Iz / h**2
K[11, 1] = 6 * E * Iz / h**2
K[5, 5] = 4 * E * Iz / h
K[11, 11] = 4 * E * Iz / h
K[5, 11] = 2 * E * Iz / h
K[11, 5] = 2 * E * Iz / h
K[5, 7] = -6 * E * Iz / h**2
K[7, 5] = -6 * E * Iz / h**2
K[11, 7] = -6 * E * Iz / h**2
K[7, 11] = -6 * E * Iz / h**2
# Bending in xz plane
K[2, 2] = 12 * E * Iy / h**3
K[8, 8] = 12 * E * Iy / h**3
K[2, 8] = -12 * E * Iy / h**3
K[8, 2] = -12 * E * Iy / h**3
K[2, 4] = -6 * E * Iy / h**2
K[4, 2] = -6 * E * Iy / h**2
K[2, 10] = -6 * E * Iy / h**2
K[10, 2] = -6 * E * Iy / h**2
K[4, 4] = 4 * E * Iy / h
K[10, 10] = 4 * E * Iy / h
K[4, 10] = 2 * E * Iy / h
K[10, 4] = 2 * E * Iy / h
K[4, 8] = 6 * E * Iy / h**2
K[8, 4] = 6 * E * Iy / h**2
K[10, 8] = 6 * E * Iy / h**2
K[8, 10] = 6 * E * Iy / h**2


i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
e_x, e_y, e_z, h = pystran.beam.beam_3d_member_geometry(i, j, xz_vector)

# Transformation matrix
T = zeros(K.shape)
T[0:3, 0] = e_x
T[0:3, 1] = e_y
T[0:3, 2] = e_z
T[3:6, 3:6] = T[0:3, 0:3]
T[6:9, 6:9] = T[0:3, 0:3]
T[9:12, 9:12] = T[0:3, 0:3]

# Transform stiffness matrix from default orientation to current orientation
K = dot(T, dot(K, T.T))

for r in range(12):
    for c in range(12):
        if abs(K[r, c] - K1[r, c]) > 1e-12 * (abs(K1[r, c]) + abs(K[r, c])):
            print(r, c, K[r, c], K1[r, c])

# plots.plot_setup(m)
# plots.plot_members(m)
# # plots.plot_deformations(m, 10.0)
# # ax = plots.plot_shear_forces(m, scale=0.50e-3)
# # ax.set_title('Shear forces')
# plots.show(m)
