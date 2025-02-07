"""
pystran - Python package for structural analysis with trusses and beams 

(C) 2025, Petr Krysl, pkrysl@ucsd.edu

# Natural Frequency of Mass supported by a Beam on Springs

Reference: Timoshenko, S., Young, D., and Weaver, W., Vibration Problems in
Engineering, John Wiley & Sons, 4th edition, 1974. page 11, problem 1.1-3.

Problem: A simple beam is supported by two spring at the endpoints. Neglecting
the distributed mass of the beam, calculate the period of free vibration of the
beam given a concentrated mass of weight W.

The answer in the book is: T = 0.533 sec., corresponding to the frequency =
1.876 CPS.
"""

from numpy.linalg import norm
from context import pystran
from pystran import model
from pystran import section
from pystran import truss
from pystran import geometry
from pystran import freedoms

# US SI(m) units
E = 2.0e5  # MPa
CTE = 1.4e-5  # 1/degC
A = 900  # mm^2
DeltaT = {1: 20.0, 2: 70.0, 3: 20.0}


def add_thermal_loads(m):
    """Set up thermal loads."""
    for member in m["truss_members"].values():
        sect = member["section"]
        EA = sect["E"] * sect["A"]
        _CTE = sect["CTE"]
        connectivity = member["connectivity"]
        _i, _j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        d = geometry.delt(_i["coordinates"], _j["coordinates"])
        nd = d / norm(d)
        _N_T = _CTE * DeltaT[member["mid"]] * EA
        model.add_load(_i, freedoms.U1, -nd[0] * _N_T)
        model.add_load(_i, freedoms.U2, -nd[1] * _N_T)
        model.add_load(_j, freedoms.U1, +nd[0] * _N_T)
        model.add_load(_j, freedoms.U2, +nd[1] * _N_T)


m = model.create(2)

model.add_joint(m, 1, [-2121.32, 2121.32])
model.add_joint(m, 2, [0.0, 2121.32])
model.add_joint(m, 3, [2121.32, 2121.32])
model.add_joint(m, 4, [0.0, 0.0])

model.add_support(m["joints"][1], freedoms.TRANSLATION_DOFS)
model.add_support(m["joints"][2], freedoms.TRANSLATION_DOFS)
model.add_support(m["joints"][3], freedoms.TRANSLATION_DOFS)


s1 = section.truss_section("s1", E=E, A=A, CTE=CTE)

model.add_truss_member(m, 1, [1, 4], s1)
model.add_truss_member(m, 2, [2, 4], s1)
model.add_truss_member(m, 3, [3, 4], s1)

# ax = plots.plot_setup(m)
# plots.plot_joint_numbers(m)
# plots.plot_members(m)
# plots.show(m)

add_thermal_loads(m)

# plots.plot_setup(m)
# plots.plot_members(m)
# ax = plots.plot_applied_forces(m, 0.001)
# ax.set_title("Forces")
# plots.show(m)

model.number_dofs(m)
model.solve_statics(m)

for j in m["joints"].values():
    print(j["jid"], ": ", j["displacements"])

for member in m["truss_members"].values():
    connectivity = member["connectivity"]
    i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
    N = truss.truss_axial_force(member, i, j)
    print("N = ", N)
    N_T = E * A * CTE * DeltaT[member["mid"]]
    print("N - N_T = ", N - N_T)


# plots.plot_setup(m)
# plots.plot_members(m)
# ax = plots.plot_deformations(m, 50.0)
# ax.set_title("Deformed shape (magnification factor = 50)")
# plots.show(m)

# plots.plot_setup(m)
# plots.plot_members(m)
# ax = plots.plot_axial_forces(m, 0.001)
# ax.set_title("Deformed shape (magnification factor = 50)")
# plots.show(m)
