"""

"""

from context import pystran
from pystran import model
from pystran import section
from pystran import plots

E = 2.0e11
G = E / (2 * (1 + 0.3))
rho = 7.85e3

h = 0.2
b = 3.0
A = b * h
Iy = b * h**3 / 12
Iz = h * b**3 / 12
Ix = Iy + Iz
J = Ix
xz_vector = [0, 0, 1]
sdeck = section.beam_3d_section(
    "sdeck", E=E, rho=rho, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)

h = 0.6
b = 0.6
A = b * h
Iy = b * h**3 / 12
Iz = h * b**3 / 12
Ix = Iy + Iz
J = Ix
xz_vector = [0, 1, 0]
spylon = section.beam_3d_section(
    "spylon", E=E, rho=rho, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)

h = 0.3
b = 0.3
A = b * h
Iy = b * h**3 / 12
Iz = h * b**3 / 12
Ix = Iy + Iz
J = Ix
xz_vector = [0, 1, 0]
scross = section.beam_3d_section(
    "scross", E=E, rho=rho, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)


h = 0.05
b = 0.05
A = b * h
scable = section.truss_section("scable", E=E, A=A, rho=rho)

m = model.create(3)

model.add_joint(m, 1, [0.0, 0.0, 0.0])
model.add_joint(m, 100, [0.0, 0.0, 1.0])
model.add_joint(m, 2, [0.0, -10.0, -0.5])
model.add_joint(m, 200, [0.0, -10.0, 0.5])
model.add_joint(m, 3, [0.0, -20.0, -1.5])
model.add_joint(m, 4, [0.0, +10.0, -0.5])
model.add_joint(m, 400, [0.0, +10.0, 0.5])
model.add_joint(m, 5, [0.0, +20.0, -1.5])

model.add_joint(m, 6, [0.0, -10.0, -4.5])
model.add_joint(m, 7, [-5.0, -10.0, 2.0])
model.add_joint(m, 8, [+5.0, -10.0, 2.0])

model.add_joint(m, 9, [0.0, +10.0, -4.5])
model.add_joint(m, 10, [-5.0, +10.0, 2.0])
model.add_joint(m, 11, [+5.0, +10.0, 2.0])

model.add_joint(m, 12, [-3.0, 0.0, 0.5])
model.add_joint(m, 13, [+3.0, 0.0, 0.5])

model.add_joint(m, 14, [-4.0, -20.0, -1.5])
model.add_joint(m, 15, [+4.0, -20.0, -1.5])

model.add_joint(m, 16, [-4.0, +20.0, -1.5])
model.add_joint(m, 17, [+4.0, +20.0, -1.5])

model.add_support(m["joints"][3], freedoms.TRANSLATION_DOFS)
model.add_support(m["joints"][5], freedoms.TRANSLATION_DOFS)
model.add_support(m["joints"][6], freedoms.ALL_DOFS)
model.add_support(m["joints"][9], freedoms.ALL_DOFS)

model.add_support(m["joints"][14], freedoms.ALL_DOFS)
model.add_support(m["joints"][15], freedoms.ALL_DOFS)
model.add_support(m["joints"][16], freedoms.ALL_DOFS)
model.add_support(m["joints"][17], freedoms.ALL_DOFS)

model.add_beam_member(m, 1, [100, 200], sdeck)
model.add_beam_member(m, 2, [200, 3], sdeck)
model.add_beam_member(m, 3, [100, 400], sdeck)
model.add_beam_member(m, 4, [400, 5], sdeck)

model.add_beam_member(m, 600, [1, 100], spylon)
model.add_beam_member(m, 601, [2, 200], spylon)
model.add_beam_member(m, 603, [4, 400], spylon)

model.add_beam_member(m, 6, [2, 6], spylon)
model.add_beam_member(m, 7, [2, 7], spylon)
model.add_beam_member(m, 8, [2, 8], spylon)
model.add_beam_member(m, 9, [4, 9], spylon)
model.add_beam_member(m, 10, [4, 10], spylon)
model.add_beam_member(m, 11, [4, 11], spylon)

model.add_beam_member(m, 12, [1, 12], scross)
model.add_beam_member(m, 13, [1, 13], scross)

model.add_truss_member(m, 1, [14, 7], scable)
model.add_truss_member(m, 2, [7, 12], scable)
model.add_truss_member(m, 3, [12, 10], scable)
model.add_truss_member(m, 4, [10, 16], scable)
model.add_truss_member(m, 5, [15, 8], scable)
model.add_truss_member(m, 6, [8, 13], scable)
model.add_truss_member(m, 7, [13, 11], scable)
model.add_truss_member(m, 8, [11, 17], scable)

model.add_load(m["joints"][1], freedoms.U3, -5e3 - 7.5e3)


model.number_dofs(m)

print("Number of free degrees of freedom = ", m["nfreedof"])
print("Number of all degrees of freedom = ", m["ntotaldof"])

# print([j["dof"] for j in m["joints"].values()])

model.solve_free_vibration(m)


# print([j["displacements"] for j in m["joints"].values()])

# print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

# print(m["U"][0 : m["nfreedof"]])


# if (
#     norm(b["displacements"] - [0.0, 0.0, -0.02238452, 0.00419677, 0.00593197, 0.0])
#     > 1.0e-5.
# ):
#     raise ValueError("Displacement calculation error")

# print("Displacement calculation OK")

# print('Reference: ', [-0.02238452,  0.00419677,  0.00593197])

plots.plot_setup(m)
plots.plot_members(m)
model.set_solution(m, m["eigvecs"][:, 0])
plots.plot_deformations(m, 500.0)
# ax = plots.plot_shear_forces(m, scale=0.50e-3)
# ax.set_title('Shear forces')
plots.show(m)
