"""
pystran - Python package for structural analysis with trusses and beams 

(C) 2025, Petr Krysl, pkrysl@ucsd.edu

# Example of a three dimensional frame problem

This example is completely solved in the book Matrix Analysis of Structures by
Robert E. Sennett, ISBN 978-1577661436 (Section 6.4). 

Displacements and internal forces are provided in the book, and we can check our
solution against these reference values.


Important note: Our orientation of the local coordinate system is such that web
of the H-beams is parallel to beam z axis! This is different from the orientation in
the book, where the web is parallel to the y axis.
"""

# We begin with the standard imports:

from math import cos, sin, pi
from numpy import array
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import geometry
from pystran import freedoms
from pystran import plots
from pystran import beam

# Important note: Our orientation of the local coordinate system is such that
# web of the H-beams is parallel to z axis! This is different from the
# orientation in the book, where the web is parallel to the y axis.

# Define a few constants:
# US customary units, inches, pounds, seconds are assumed.
E = 29e6

kbars = 100000.0  # lb/in

ksprings = kbars * 1000000.0  # lb/in


def compute_area(_m, connectivity, _E, _kbars):
    _i, _j = _m["joints"][connectivity[0]], _m["joints"][connectivity[1]]
    L = geometry.vlen(_i["coordinates"], _j["coordinates"])
    return _kbars * L / _E


# The model is created as three dimensional (3 as argument!).
m = model.create(2)

# Joints are added at their locations, and the supports of the clamped joints are specified.
model.add_joint(m, 1, [-20 * 12, -10 * 12])
model.add_joint(m, 2, [-10 * 12, 0.0])
model.add_joint(m, 3, [0.0, -10 * 12])


# There are three truss bars. The cross sectional properties are compute so that
# the axial stiffness is the same for all bars.

connectivity = [1, 2]
model.add_truss_member(
    m,
    1,
    connectivity,
    section.truss_section("s1", E=E, A=compute_area(m, connectivity, E, kbars)),
)
connectivity = [1, 3]
model.add_truss_member(
    m,
    2,
    connectivity,
    section.truss_section("s2", E=E, A=compute_area(m, connectivity, E, kbars)),
)
connectivity = [3, 2]
model.add_truss_member(
    m,
    3,
    connectivity,
    section.truss_section("s3", E=E, A=compute_area(m, connectivity, E, kbars)),
)
print(m)


# Next we add the applied moment, and
model.add_load(m["joints"][2], freedoms.U1, 10)
model.add_load(m["joints"][2], freedoms.U2, 5)


model.add_extension_spring_to_ground(m["joints"][1], 1, [0, 1], ksprings)
model.add_extension_spring_to_ground(m["joints"][1], 2, [1, 0], ksprings)
model.add_extension_spring_to_ground(
    m["joints"][3], 1, [-cos(30 / 180 * pi), sin(30 / 180 * pi)], ksprings
)

# Now we can solve the static equilibrium of the frame.
model.number_dofs(m)
model.solve_statics(m)

# The displacements of the joints can be printed out. Recall that joints 3 and 4
# are clamped, and their displacements are therefore zero.
for jid in [1, 2, 3]:
    j = m["joints"][jid]
    print(jid, kbars * j["displacements"])

ax = plots.plot_setup(m)
plots.plot_members(m)
ax = plots.plot_deformations(m, 100000.0)
ax.set_title("Deformations (x100000)")
plots.show(m)
