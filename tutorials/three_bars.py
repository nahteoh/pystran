"""
Analysis of Geometrically 
Nonlinear Structures 
Second Edition 

by 
Robert Levy 
Technion-Israel Institute of Technology, 
Haifa, Israel 

and 
William R. Spillers 
New Jersey Institute of Technology, 

Three bar example on page 32
"""

# Due to the Python source code loading conventions, pystran is not immediately
# visible to the Python interpreter. Loading it through the context is a
# workaround.
from context import pystran

# The following imports are necessary for the example to work.
from numpy.linalg import norm
from numpy import concatenate, dot
from pystran import model
from pystran import section
from pystran import truss
from pystran import plots

# Now we create the model `m`, which is a dictionary that contains all the
# constituent parts of the system. The model is created for a planar structure
# (2). The deformation is assumed to occur in the plane `x-z`.
m = model.create(2)

# There are four joints. Note that the coordinates are supplied without physical
# units. Everything needs to be provided in consistent units.
model.add_joint(m, 1, [10.0, 20.0])
model.add_joint(m, 2, [0.0, 20.0])
model.add_joint(m, 3, [0.0, 10.0])
model.add_joint(m, 4, [+10.0, 0.0])

# For instance, the material parameters and the cross sectional area of the bars.
E = 30000000.0
A = 0.65700000

# These are used to create a section.
s1 = section.truss_section("steel", E, A)

# The truss members are added to the model. The connectivity is provided as
# lists of [start, end] joints. The section is also provided.
model.add_truss_member(m, 1, [1, 2], s1)
model.add_truss_member(m, 2, [1, 3], s1)
model.add_truss_member(m, 3, [1, 4], s1)

# The supports are added to the model. The pinned supports are added to the
# joints 2, 3, and 4, by constraining both x and z displacements to zero.
for i in [2, 3, 4]:
    for d in range(2):
        model.add_support(m["joints"][i], d)

# Loads are added to the joint 1. The loads are provided as [x, z] components.
model.add_load(m["joints"][1], 0, -10000.0)
model.add_load(m["joints"][1], 1, -10000.0 / 2.0)

# The model is shown graphically, members are displayed with the member numbers attached.
plots.plot_setup(m)
plots.plot_members(m)
plots.plot_member_numbers(m)
plots.show(m)

# The degrees of freedom are numbered, first the unknowns (free) degrees of
# freedom, then the known (prescribed) degrees of freedom.
model.number_dofs(m)
print("Total Degrees of Freedom = ", m["ntotaldof"])
print("Free Degrees of Freedom = ", m["nfreedof"])

# The model is solved.  This means that the displacements are calculated  from
# the balance of the joints.
model.solve(m)

# The displacements are printed.
for j in m["joints"].values():
    print(j["displacements"])

# The displacements are compared to the reference values from the book.
if norm(m["joints"][1]["displacements"] - [-0.0033325938, -0.001591621]) > 1.0e-3:
    raise ValueError("Displacement calculation error")
print("Displacement calculation OK")

# The forces in the bars can be calculated using the strain-displacement matrix
# to compute the strain, from which we can evaluate the force.
for b in m["truss_members"].values():
    connectivity = b["connectivity"]
    i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
    e_x, L = truss.truss_member_geometry(i, j)
    B = truss.truss_strain_displacement(e_x, L)
    u = concatenate((i["displacements"], j["displacements"]))
    eps = dot(B, u)
    print("Bar " + str(connectivity) + " force = ", E * A * eps[0])

# The forces in the bars can be calculated using the `truss_axial_force` function.
# That function simply does what was described in the loop above.
for b in m["truss_members"].values():
    connectivity = b["connectivity"]
    i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
    N = truss.truss_axial_force(b, i, j)
    print("N = ", N)

# These are the reference values of the forces from the book.
print("Reference forces: ", -0.656854250e4, -0.48528137e4, -0.15685425e4)

# The solution is visualized with deformed shape.
plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_deformations(m, 100.0)
ax.set_title("Deformed shape (magnification factor = 100)")
plots.show(m)


# The solution is further visualized with graphical representation of the
# internal (axial) forces.
plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_axial_forces(m, 1 / 3000.0)
ax.set_title("Axial Forces")
plots.show(m)
