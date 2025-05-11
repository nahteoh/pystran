from numpy import dot
from numpy.linalg import norm
import context
from pystran import model
from pystran import section
from pystran import freedoms
from pystran import plots



# Create a two dimensional (planar) model.
m = model.create(2)

# There are three joints. (in m)
model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [6, 0.0])
model.add_joint(m, 3, [12, 0.0])
model.add_joint(m, 4, [18, 0.0])


model.add_support(m["joints"][1], freedoms.U1)
model.add_support(m["joints"][1], freedoms.U2)
model.add_support(m["joints"][1], freedoms.UR3)
model.add_support(m["joints"][2], freedoms.U2)
model.add_support(m["joints"][3], freedoms.U2)
model.add_support(m["joints"][4], freedoms.U1)
model.add_support(m["joints"][4], freedoms.U2)
model.add_support(m["joints"][4], freedoms.UR3)



E = 70000000 # kn/m^2 
I = 0.000225 # m^4
A =  0.0476222# m^2
s1 = section.beam_2d_section("section_1", E, A , I)

model.add_beam_member(m, 1, [1, 2], s1)
model.add_beam_member(m, 2, [2, 3], s1)
model.add_beam_member(m, 3, [3, 4], s1)


model.add_load(m["joints"][2], freedoms.UR3, -150)

model.number_dofs(m)
nt = m["ntotaldof"]
nf = m["nfreedof"]
print("Total Degrees of Freedom = ", nt)
print("Free Degrees of Freedom = ", nf)



model.solve_statics(m)
print(m["joints"])
# The stiffness matrix for the free degrees of freedom can be printed.
print(m["K"][0:3, 0:3])

# Here are the calculated free degrees of freedom:
print(m["U"][0:3])




# We can calculate the reactions at the supports. This is the manual approach
# to that using the partitioning of the stiffness matrix and the displacement
# vector.
Kdf = m["K"][nf:nt, 0:nf]
Kdd = m["K"][nf:nt, nf:nt]
Uf = m["U"][0:nf]
Ud = m["U"][nf:nt]

print('Displacement:', Uf)
# The reactions follow:
R = dot(Kdf, Uf) + dot(Kdd, Ud)
print("Reactions = ", R)

# In order to understand moment and shear diagrams, we start by plotting the
# geometry with the orientation of the local coordinate system on each beam.
plots.setup(m)
plots.plot_members(m)
ax = plots.plot_member_orientation(m, 1.0)
ax.set_title("Local coordinate systems (red -- local x, blue -- local z)")
plots.show(m)


# The deformed shape shows the curvatures of the beam.
plots.setup(m)
plots.plot_members(m)
ax = plots.plot_deformations(m, 100.0)
ax.set_title("Deformed shape (magnified 100 times)")
plots.show(m)

# The bending moment can be compared with the curvature of the beam.
plots.setup(m)
plots.plot_members(m)
ax = plots.plot_bending_moments(m, scale=0.1)
ax.set_title("Moments")
plots.show(m)

# The shear forces are the slopes of the moment diagram.
plots.setup(m)
plots.plot_members(m)
ax = plots.plot_shear_forces(m)
ax.set_title("Shear forces")
plots.show(m)

# The reaction forces can be plotted after the reactions have been computed.
model.statics_reactions(m)
plots.setup(m)
plots.plot_members(m)
ax = plots.plot_reaction_forces(m)
ax.set_title("Reaction forces")
plots.show(m)
