"""
Simple cantilever under a concentrated force at the end.
Calculate the stiffness matrix.
"""

from context import pystran
from pystran import model
from pystran import section
from pystran import geometry
from pystran import plots
from math import sqrt

m = model.create(2)

h = 8.0
model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [h, 0.0])

model.add_support(m["joints"][1], freedoms.ALL_DOFS)
model.add_support(m["joints"][2], freedoms.U1)


E = 2.0e11
A = 6000 / 10**6
I = 200e6 / 10**12
s1 = section.beam_2d_section("sect_1", E=E, A=A, I=I)
model.add_beam_member(m, 1, [1, 2], s1)

F = 10e3
model.add_load(m["joints"][2], freedoms.U2, F)

model.number_dofs(m)

model.solve_statics(m)

print(f"Stiffness = {m['K']}")
nf = m["nfreedof"]
Kff = m["K"][0:nf, 0:nf]
print(f"Kff = {Kff}")
print(f"12*E*I/h**3 = {12*E*I/h**3}")
print(f"6*E*I/h**2 = {6*E*I/h**2}")
print(f"4*E*I/h = {4*E*I/h}")
