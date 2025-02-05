"""

"""

from context import pystran
from pystran import model
from pystran import section
from pystran import geometry
from pystran import plots
from numpy import dot
from numpy.linalg import norm

E = 210e9
nu = 0.3
G = E / (2 * (1 + nu))

A, Ix, Iy, Iz, J = section.square_tube(0.1, 0.1, 0.008, 0.008)

m = model.create(2)

model.add_joint(m, 1, [-0.5, 0.0])
model.add_joint(m, 2, [0.0, 0.0])
model.add_joint(m, 3, [0.5, 0.0])
model.add_support(m["joints"][1], freedoms.ALL_DOFS)
model.add_support(m["joints"][3], freedoms.ALL_DOFS)

s1 = section.beam_2d_section("s1", E, A, Iy)
model.add_beam_member(m, 1, [1, 2], s1)
model.add_beam_member(m, 2, [2, 3], s1)


model.add_load(m["joints"][2], freedoms.U2, -100e3)

ax = plots.plot_setup(m)
plots.plot_joint_numbers(m)
plots.plot_translation_supports(m, scale=0.001)
ax.set_title("Translation supports")
plots.show(m)


ax = plots.plot_setup(m)
plots.plot_joint_numbers(m)
plots.plot_rotation_supports(m, scale=0.001)
ax.set_title("Rotation supports")
plots.show(m)


model.number_dofs(m)

model.solve_statics(m)

print([j["displacements"] for j in m["joints"].values()])

# if norm(m['U'][0:3] - [-0.02969075, -0.02742406, 0.03952194]) > 1.e-3:
#     raise ValueError('Displacement calculation error')
# else:
#     print('Displacement calculation OK')


plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 1000.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_bending_moments(m, 0.0001)
ax.set_title("Moments")
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_shear_forces(m, 0.0001)
ax.set_title("Shear forces")
plots.show(m)


m = model.create(2)

model.add_joint(m, 1, [-0.75, 0.0])
model.add_joint(m, 2, [0.0, 0.0])
model.add_joint(m, 3, [0.25, 0.0])
model.add_support(m["joints"][1], freedoms.ALL_DOFS)
model.add_support(m["joints"][3], freedoms.ALL_DOFS)

s1 = section.beam_2d_section("s1", E, A, Iy)
model.add_beam_member(m, 1, [1, 2], s1)
model.add_beam_member(m, 2, [2, 3], s1)


model.add_load(m["joints"][2], freedoms.U2, -100e3)

# ax = plots.plot_setup(m)
# plots.plot_joint_numbers(m)
# plots.plot_translation_supports(m, scale=0.001)
# ax.set_title("Translation supports")
# plots.show(m)

model.number_dofs(m)

model.solve_statics(m)

print([j["displacements"] for j in m["joints"].values()])

# if norm(m['U'][0:3] - [-0.02969075, -0.02742406, 0.03952194]) > 1.e-3:
#     raise ValueError('Displacement calculation error')
# else:
#     print('Displacement calculation OK')

# The reactions at the supports can be calculated. The reactions at joints 3 and 4
# are compared to the reference values.
model.statics_reactions(m)

for jid in [1, 3]:
    j = m["joints"][jid]
    print(f"Joint {jid}:")
    print(
        f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, My={j['reactions'][2]:.5}"
    )

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 1000.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_bending_moments(m, 0.0001)
ax.set_title("Moments")
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_shear_forces(m, 0.0001)
ax.set_title("Shear forces")
plots.show(m)


m = model.create(3)

model.add_joint(m, 1, [-0.75, 0.0, 0.0])
model.add_joint(m, 2, [0.0, 0.0, 0.0])
model.add_joint(m, 3, [0.25, 0.0, 0.0])
model.add_support(m["joints"][1], freedoms.ALL_DOFS)
model.add_support(m["joints"][3], freedoms.ALL_DOFS)

xz_vector = [0, 0, 1]
s1 = section.beam_3d_section(
    "s1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
)
model.add_beam_member(m, 1, [1, 2], s1)
model.add_beam_member(m, 2, [2, 3], s1)


model.add_load(m["joints"][2], freedoms.U2, -100e3)

ax = plots.plot_setup(m)
plots.plot_joint_numbers(m)
plots.plot_translation_supports(m, scale=0.001)
ax.set_title("Translation supports")
plots.show(m)

ax = plots.plot_setup(m)
plots.plot_joint_numbers(m)
plots.plot_rotation_supports(m, scale=0.001)
ax.set_title("Rotation supports")
plots.show(m)

model.number_dofs(m)

model.solve_statics(m)

print([j["displacements"] for j in m["joints"].values()])

# if norm(m['U'][0:3] - [-0.02969075, -0.02742406, 0.03952194]) > 1.e-3:
#     raise ValueError('Displacement calculation error')
# else:
#     print('Displacement calculation OK')

# The reactions at the supports can be calculated. The reactions at joints 3 and 4
# are compared to the reference values.
model.statics_reactions(m)

for jid in [1, 3]:
    j = m["joints"][jid]
    print(f"Joint {jid}:")
    print(
        f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, My={j['reactions'][2]:.5}"
    )

plots.plot_setup(m)
plots.plot_members(m)
plots.plot_deformations(m, 1000.0)
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_bending_moments(m, 0.0001, "z")
ax.set_title("Moments")
plots.show(m)

plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_shear_forces(m, 0.0001, "y")
ax.set_title("Shear forces")
plots.show(m)
