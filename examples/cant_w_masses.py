"""
Created on 01/19/2025

Structural Analysis: A Unified Classical and Matrix, Ghali, Amin; Neville, Adam
-- Edition 7, 2017, Taylor and Francis

Example 24.2 - Cantilever with added masses

The mass density of the beam itself is artificially reduced so that there are
only the added masses.
"""

from math import sqrt, pi
from numpy import array
from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import freedoms
from pystran import plots

E = 2.0e11
G = E / (2 * (1 + 0.3))
rho = 7.85e3 / 10000  # artificially reduce the mass density of the beam

h = 0.12
b = 0.03
A = b * h
Iy = b * h**3 / 12
sbar = section.beam_2d_section("sbar", E=E, rho=rho, A=A, I=Iy)
L = 2.0
W = 3.0 * 9.81
g = 9.81


m = model.create(2)

model.add_joint(m, 1, [0.0, 3 * L])
model.add_joint(m, 2, [0.0, 2 * L])
model.add_joint(m, 3, [0.0, 1 * L])
model.add_joint(m, 4, [0.0, 0.0])

model.add_support(m["joints"][4], freedoms.ALL_DOFS)

model.add_beam_member(m, 1, [1, 2], sbar)
model.add_beam_member(m, 2, [2, 3], sbar)
model.add_beam_member(m, 3, [3, 4], sbar)

model.add_mass(m["joints"][1], freedoms.U1, 4 * W / g)
model.add_mass(m["joints"][1], freedoms.U2, 4 * W / g)
model.add_mass(m["joints"][2], freedoms.U1, W / g)
model.add_mass(m["joints"][2], freedoms.U2, W / g)
model.add_mass(m["joints"][3], freedoms.U1, W / g)
model.add_mass(m["joints"][3], freedoms.U2, W / g)

model.number_dofs(m)

model.solve_free_vibration(m)

expected = array([0.1609, 1.7604, 5.0886]) * sqrt(g * E * Iy / W / L**3) / 2 / pi
print("Expected frequencies (zero mass of beam): ", expected)
print("Computed frequencies: ", m["frequencies"][0:3])
if norm((m["frequencies"][0:3] - expected) / expected) > 1.0e-2:
    raise ValueError("Frequency calculation error")

for mode in range(3):
    plots.plot_setup(m)
    plots.plot_members(m)
    model.set_solution(m, m["eigvecs"][:, mode])
    ax = plots.plot_deformations(m, 50.0)
    ax.set_title(f"Mode {mode}: f = {sqrt(m['eigvals'][mode])/2/pi:.3f} Hz")
    plots.show(m)
