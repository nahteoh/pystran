"""
Created on 01/12/2025

Example 4.13 from 
Matrix Structural Analysis: Second Edition 2nd Edition
by William McGuire, Richard H. Gallagher, Ronald D. Ziemian 
"""

from numpy import array, dot
from numpy.linalg import cross
from context import pystran
from pystran import model
from pystran import section

h = 8.0
E = 2.0e6
G = E / (2 * (1 + 0.3))
H = 0.13
B = 0.5
A = H * B
Iy = H * B**3 / 12
Iz = H**3 * B / 12
Ix = Iy + Iz
J = Ix


def test(e_x, e_y, e_z, F, refdefl, refslope):
    m = model.create(3)
    model.add_joint(m, 1, [0.0, 0.0, 0.0])
    model.add_joint(m, 2, h * e_x)
    clamped = m["joints"][1]
    freeend = m["joints"][2]

    model.add_support(clamped, model.U1)
    model.add_support(clamped, model.U2)
    model.add_support(clamped, model.U3)
    model.add_support(clamped, model.UR1)
    model.add_support(clamped, model.UR2)
    model.add_support(clamped, model.UR3)

    s1 = section.beam_3d_section(
        "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=e_z
    )
    model.add_beam_member(m, 1, [1, 2], s1)

    model.add_load(freeend, model.U1, F[0])
    model.add_load(freeend, model.U2, F[1])
    model.add_load(freeend, model.U3, F[2])

    model.number_dofs(m)

    # print('Number of free degrees of freedom = ', m['nfreedof'])
    # print('Number of all degrees of freedom = ', m['ntotaldof'])

    # print([j['dof'] for j in m['joints'].values()])

    model.solve(m)

    # print([j['displacements'] for j in m['joints'].values()])

    # print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

    # print(m['U'][0:m['nfreedof']])

    rot = cross(e_x, F)
    defl = dot(F, freeend["displacements"][0:3])
    slope = dot(rot, freeend["displacements"][3:6])

    print("Reference deflection: ", refdefl, " Computed deflection: ", defl)
    print("Reference slope: ", refslope, " Computed slope: ", slope)

    if abs(defl - refdefl) > 1.0e-3 * abs(refdefl):
        raise ValueError("Displacement calculation error")

    if abs(slope - refslope) > 1.0e-3 * abs(refslope):
        raise ValueError("Slope calculation error")

    # plots.plot_setup(m)
    # plots.plot_members(m)
    # plots.plot_deformations(m, 10.0)
    # # ax = plots.plot_shear_forces(m, scale=0.50e-3)
    # # ax.set_title('Shear forces')
    # plots.show(m)


e_x, e_y, e_z = array([1, 0, 0]), array([0, 1, 0]), array([0, 0, 1])
F = array([0, 1, 0])
refdefl = 1 * h**3 / (3 * E * Iz)
refslope = 1 * h**2 / (2 * E * Iz)
test(e_x, e_y, e_z, F, refdefl, refslope)

e_x, e_y, e_z = array([1, 0, 0]), array([0, 1, 0]), array([0, 0, 1])
F = array([0, 0, 1])
refdefl = 1 * h**3 / (3 * E * Iy)
refslope = 1 * h**2 / (2 * E * Iy)
test(e_x, e_y, e_z, F, refdefl, refslope)

e_x, e_y, e_z = array([0, 1, 0]), array([0, 0, 1]), array([1, 0, 0])
F = array([0, 0, 1])
refdefl = 1 * h**3 / (3 * E * Iz)
refslope = 1 * h**2 / (2 * E * Iz)
test(e_x, e_y, e_z, F, refdefl, refslope)

e_x, e_y, e_z = array([0, 1, 0]), array([0, 0, 1]), array([1, 0, 0])
F = array([0, 0, -1])
refdefl = 1 * h**3 / (3 * E * Iz)
refslope = 1 * h**2 / (2 * E * Iz)
test(e_x, e_y, e_z, F, refdefl, refslope)
