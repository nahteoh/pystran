"""
pystran unit tests
"""

import unittest

from numpy import array, dot, outer, concatenate
from numpy.linalg import norm
from pystran import model
from pystran import section
from pystran import plots
from pystran import beam
from pystran import rotation


class UnitTestsPlaneTrusses(unittest.TestCase):

    def test_truss_dome(self):
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


        Space truss dome in section 2.4.2

        Vertical deflection at the crown: -.20641184e+00 in (linear analysis)
        """

        m = model.create(3)

        model.add_joint(m, 1, [0.0, 0.0, 0.32346000e1])
        model.add_joint(m, 2, [0.49212500e1, 0.85239000e1, 0.24472000e1])
        model.add_joint(m, 3, [-0.49212500e1, 0.85239000e1, 0.24472000e1])
        model.add_joint(m, 4, [-0.98425000e1, 0.0, 0.24472000e1])
        model.add_joint(m, 5, [-0.49212500e1, -0.85239000e1, 0.24472000e1])
        model.add_joint(m, 6, [0.49212500e1, -0.85239000e1, 0.24472000e1])
        model.add_joint(m, 7, [0.98425000e1, 0.0, 0.24472000e1])
        model.add_joint(m, 8, [0.0, 0.19685000e02, 0.0])
        model.add_joint(m, 9, [-0.17047200e02, 0.98425000e1, 0.0])
        model.add_joint(m, 10, [-0.17047200e02, -0.98425000e1, 0.0])
        model.add_joint(m, 11, [0.0, -0.19685000e02, 0.0])
        model.add_joint(m, 12, [0.17047200e02, -0.98425000e1, 0.0])
        model.add_joint(m, 13, [0.17047200e02, 0.98425000e1, 0.0])

        E = 30000000.0
        A = 0.0155
        s1 = section.truss_section("steel", E, A)

        for id, j in enumerate(
            [
                [1, 2],
                [1, 3],
                [1, 4],
                [1, 5],
                [1, 6],
                [1, 7],
                [6, 12],
                [7, 12],
                [7, 13],
                [2, 13],
                [2, 8],
                [3, 8],
                [3, 9],
                [4, 9],
                [4, 10],
                [5, 10],
                [5, 11],
                [6, 11],
                [2, 3],
                [3, 4],
                [4, 5],
                [5, 6],
                [6, 7],
                [7, 2],
            ]
        ):
            model.add_truss_member(m, id, j, s1)

        for i in [8, 9, 10, 11, 12, 13]:
            for d in range(m["dim"]):
                model.add_support(m["joints"][i], d)

        model.add_load(m["joints"][1], 2, -220.46)

        model.number_dofs(m)
        print("Total Degrees of Freedom = ", m["ntotaldof"])
        print("Free Degrees of Freedom = ", m["nfreedof"])

        model.solve(m)

        for j in m["joints"].values():
            print(j["displacements"])

        if norm(m["joints"][1]["displacements"][2] - (-0.20641184e00)) > 1.0e-3 * abs(
            -0.20641184e00
        ):
            raise ValueError("Displacement calculation error")
        else:
            print("Displacement calculation OK")

        # for b in m['truss_members'].values():
        #     connectivity = b['connectivity']
        #     i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        #     e_x, L = truss.truss_member_geometry(i, j)
        #     B = truss.strain_displacement(e_x, L)
        #     u = concatenate((i['displacements'], j['displacements']))
        #     eps = dot(B, u)
        #     print('Bar ' + str(connectivity) + ' force = ', E * A * eps[0])

        # plots.plot_setup(m)
        # plots.plot_members(m)
        # plots.plot_deformations(m, 5.0)
        # plots.show(m)


class UnitTestsSpaceFrames(unittest.TestCase):

    def test_spfr_weaver_gere_1(self):
        """
        Created on 01/12/2025

        Example 5.8 from Matrix Structural Analysis: Second Edition 2nd Edition by
        William McGuire, Richard H. Gallagher, Ronald D. Ziemian

        The section properties are not completely defined in the book.  They are
        taken from example 4.8, which does not provide both second moments of area.
        They are taken here as both the same.
        """

        # SI units
        L = 3.0
        E = 200e9
        G = 80e9
        A = 0.01
        Iz = 1e-3
        Iy = 1e-3
        Ix = 2e-3
        J = Ix  # Torsional constant
        P = 60000

        m = model.create(3)

        jA, jB, jC, jD, jE = 3, 1, 2, 4, 5
        model.add_joint(m, jA, [0.0, 0.0, 0.0])
        model.add_joint(m, jB, [0.0, L, 0.0])
        model.add_joint(m, jC, [2 * L, L, 0.0])
        model.add_joint(m, jD, [3 * L, 0.0, L])
        model.add_joint(m, jE, [L, L, 0.0])

        model.add_support(m["joints"][jA], model.CLAMPED)
        model.add_support(m["joints"][jD], model.CLAMPED)

        xz_vector = [1, 0, 0]
        s1 = section.beam_3d_section(
            "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )
        model.add_beam_member(m, 1, [jA, jB], s1)
        xz_vector = [0, 1, 0]
        s2 = section.beam_3d_section(
            "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )
        model.add_beam_member(m, 2, [jE, jB], s2)
        model.add_beam_member(m, 3, [jE, jC], s2)
        model.add_beam_member(m, 4, [jC, jD], s1)

        model.add_load(m["joints"][jB], model.U1, 2 * P)
        model.add_load(m["joints"][jE], model.U3, 4 * P)
        model.add_load(m["joints"][jC], model.U2, -P)
        model.add_load(m["joints"][jC], model.UR3, -P * L)

        model.number_dofs(m)

        # print("Number of free degrees of freedom = ", m["nfreedof"])
        # print("Number of all degrees of freedom = ", m["ntotaldof"])

        # print([j['dof'] for j in m['joints'].values()])

        model.solve(m)

        for jid in [jB, jC, jE]:
            j = m["joints"][jid]
            print(jid, j["displacements"])

        # print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

        # print(m['U'][0:m['nfreedof']])

        if (
            norm(
                m["joints"][1]["displacements"]
                - [
                    -8.59409726e-04,
                    5.77635277e-05,
                    5.00764459e-03,
                    2.39333188e-03,
                    -1.62316861e-03,
                    6.81331291e-04,
                ]
            )
            > 1.0e-5
        ):
            raise ValueError("Displacement calculation error")
        else:
            print("Displacement calculation OK")

        if (
            norm(
                m["joints"][2]["displacements"]
                - [
                    -0.00117605,
                    0.00325316,
                    0.00525552,
                    0.00128843,
                    0.00172094,
                    -0.00077147,
                ]
            )
            > 1.0e-5
        ):
            raise ValueError("Displacement calculation error")
        else:
            print("Displacement calculation OK")

        # print('Reference: ', [-0.02238452,  0.00419677,  0.00593197])

        # plots.plot_setup(m)
        # plots.plot_members(m)
        # # plots.plot_member_numbers(m)
        # plots.plot_deformations(m, 100.0)
        # # ax = plots.plot_shear_forces(m, scale=0.50e-3)
        # # ax.set_title('Shear forces')
        # plots.show(m)

    def test_weaver_1(self):
        """
        Created on 01/20/2025

        Weaver Jr., W., Computer Programs for Structural Analysis, page 146,
        problem 8. From: STAAD.Pro 2023.00.03
        """

        # US customary units, inches, pounds, seconds
        L = 120.0
        E = 30000
        G = E / (2 * (1 + 0.3))
        A = 11
        Iz = 56
        Iy = 56
        Ix = 83
        J = Ix  # Torsional constant
        F = 2
        P = 1
        M = 120

        m = model.create(3)

        model.add_joint(m, 3, [0.0, 0.0, 0.0])
        model.add_joint(m, 1, [0.0, L, 0.0])
        model.add_joint(m, 2, [2 * L, L, 0.0])
        model.add_joint(m, 4, [3 * L, 0.0, L])

        model.add_support(m["joints"][3], model.CLAMPED)
        model.add_support(m["joints"][4], model.CLAMPED)

        xz_vector = [1, 0, 0]
        sect_1 = section.beam_3d_section(
            "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )
        model.add_beam_member(m, 1, [3, 1], sect_1)
        xz_vector = [0, 1, 0]
        sect_2 = section.beam_3d_section(
            "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )
        model.add_beam_member(m, 2, [1, 2], sect_2)
        model.add_beam_member(m, 3, [2, 4], sect_2)

        model.add_load(m["joints"][1], model.U1, F)
        model.add_load(m["joints"][2], model.U2, -P)
        model.add_load(m["joints"][2], model.UR3, -M)

        model.number_dofs(m)

        # print("Number of free degrees of freedom = ", m["nfreedof"])
        # print("Number of all degrees of freedom = ", m["ntotaldof"])

        # print([j['dof'] for j in m['joints'].values()])

        model.solve(m)

        for jid in [1, 2]:
            j = m["joints"][jid]
            print(jid, j["displacements"])

        # print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

        # print(m['U'][0:m['nfreedof']])

        # 0.22267   0.00016  -0.17182  -0.00255   0.00217  -0.00213
        # 0.22202  -0.48119  -0.70161  -0.00802   0.00101  -0.00435
        ref1 = [0.22267, 0.00016, -0.17182, -0.00255, 0.00217, -0.00213]
        if norm(m["joints"][1]["displacements"] - ref1) > 1.0e-1 * norm(ref1):
            raise ValueError("Displacement calculation error")
        else:
            print("Displacement calculation OK")
        ref2 = [0.22202, -0.48119, -0.70161, -0.00802, 0.00101, -0.00435]
        if norm(m["joints"][2]["displacements"] - ref2) > 1.0e-1 * norm(ref2):
            raise ValueError("Displacement calculation error")
        else:
            print("Displacement calculation OK")

        # plots.plot_setup(m)
        # plots.plot_members(m)
        # # plots.plot_member_numbers(m)
        # plots.plot_deformations(m, 30.0)
        # # ax = plots.plot_shear_forces(m, scale=0.50e-3)
        # # ax.set_title('Shear forces')
        # plots.show(m)

    def test_original_weaver_1(self):
        """
        Created on 01/12/2025

        Example 5.8 from Matrix Structural Analysis: Second Edition 2nd Edition by
        William McGuire, Richard H. Gallagher, Ronald D. Ziemian

        The section properties are not completely defined in the book.  They are
        taken from example 4.8, which does not provide both second moments of area.
        They are taken here as both the same.
        """

        # US customary units, inches, pounds, seconds
        L = 120.0
        E = 30000
        G = E / (2 * (1 + 0.3))
        A = 11
        Iz = 56
        Iy = 56
        Ix = 83
        J = Ix  # Torsional constant
        F = 2
        P = 1
        M = 120

        m = model.create(3)

        model.add_joint(m, 3, [0.0, 0.0, 0.0])
        model.add_joint(m, 1, [0.0, L, 0.0])
        model.add_joint(m, 2, [2 * L, L, 0.0])
        model.add_joint(m, 4, [3 * L, 0.0, L])

        model.add_support(m["joints"][3], model.CLAMPED)
        model.add_support(m["joints"][4], model.CLAMPED)

        xz_vector = [0, 0, 1]
        sect_1 = section.beam_3d_section(
            "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )
        xz_vector = [0, 0, 1]
        sect_2 = section.beam_3d_section(
            "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )
        xz_vector = rotation.rotate(m["joints"][2], m["joints"][4], [0, 1, 0], 90)
        sect_3 = section.beam_3d_section(
            "sect_3", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )

        model.add_beam_member(m, 1, [1, 2], sect_1)
        model.add_beam_member(m, 2, [3, 1], sect_2)
        model.add_beam_member(m, 3, [2, 4], sect_3)

        model.add_load(m["joints"][1], model.U1, F)
        model.add_load(m["joints"][2], model.U2, -P)
        model.add_load(m["joints"][2], model.UR3, -M)

        model.number_dofs(m)

        # print("Number of free degrees of freedom = ", m["nfreedof"])
        # print("Number of all degrees of freedom = ", m["ntotaldof"])

        # print([j['dof'] for j in m['joints'].values()])

        model.solve(m)

        for jid in [1, 2]:
            j = m["joints"][jid]
            print(jid, j["displacements"])

        # print(m['K'][0:m['nfreedof'], 0:m['nfreedof']])

        # print(m['U'][0:m['nfreedof']])

        # 0.22267   0.00016  -0.17182  -0.00255   0.00217  -0.00213
        # 0.22202  -0.48119  -0.70161  -0.00802   0.00101  -0.00435
        ref1 = [0.22267, 0.00016, -0.17182, -0.00255, 0.00217, -0.00213]
        if norm(m["joints"][1]["displacements"] - ref1) > 1.0e-1 * norm(ref1):
            raise ValueError("Displacement calculation error")
        else:
            print("Displacement calculation OK")
        ref2 = [0.22202, -0.48119, -0.70161, -0.00802, 0.00101, -0.00435]
        if norm(m["joints"][2]["displacements"] - ref2) > 1.0e-1 * norm(ref2):
            raise ValueError("Displacement calculation error")
        else:
            print("Displacement calculation OK")

        # for k in m["beam_members"].keys():
        #     member = m["beam_members"][k]
        #     connectivity = member["connectivity"]
        #     i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        #     f = beam.beam_3d_end_forces(member, i, j)
        #     print(f"Member {k}: ")
        #     print(
        #         f"   Joint {connectivity[0]}: N={f['Ni']:.5}, Qy={f['Qyi']:.5}, Qz={f['Qzi']:.5}, T={f['Ti']:.5}, My={f['Myi']:.5}, Mz={f['Mzi']:.5}: "
        #     )
        #     print(
        #         f"   Joint {connectivity[1]}: N={f['Nj']:.5}, Qy={f['Qyj']:.5}, Qz={f['Qzj']:.5}, T={f['Tj']:.5}, My={f['Myj']:.5}, Mz={f['Mzj']:.5}: "
        #     )
        member = m["beam_members"][1]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        f = beam.beam_3d_end_forces(member, i, j)
        # Member 1:
        self.assertAlmostEqual(f["Ni"], -0.89607, places=3)
        self.assertAlmostEqual(f["Qyi"], 0.43381, places=3)
        self.assertAlmostEqual(f["Qzi"], -0.21393, places=3)
        self.assertAlmostEqual(f["Ti"], -22.661, places=3)
        self.assertAlmostEqual(f["Myi"], 17.656, places=3)
        self.assertAlmostEqual(f["Mzi"], 36.188, places=3)
        # Nj=0.89607, Qy=-0.43381, Qz=0.21393, T=22.661, My=33.689, Mz=67.927
        # Member 2:
        member = m["beam_members"][2]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        f = beam.beam_3d_end_forces(member, i, j)
        self.assertAlmostEqual(f["Ni"], 0.43381, places=3)
        self.assertAlmostEqual(f["Qyi"], -1.1039, places=3)
        self.assertAlmostEqual(f["Qzi"], -0.21393, places=3)
        self.assertAlmostEqual(f["Ti"], 17.656, places=3)
        self.assertAlmostEqual(f["Myi"], 48.333, places=3)
        self.assertAlmostEqual(f["Mzi"], -96.284, places=3)
        # Nj=-0.43381, Qy=1.1039, Qz=0.21393, T=-17.656, My=-22.661, Mz=-36.188
        # Member 3:
        member = m["beam_members"][3]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        f = beam.beam_3d_end_forces(member, i, j)
        self.assertAlmostEqual(f["Ni"], -1.4687, places=3)
        self.assertAlmostEqual(f["Qyi"], 0.71755, places=3)
        self.assertAlmostEqual(f["Qzi"], 0.48234, places=3)
        self.assertAlmostEqual(f["Ti"], 36.432, places=3)
        self.assertAlmostEqual(f["Myi"], -15.499, places=3)
        self.assertAlmostEqual(f["Mzi"], 52.845, places=3)
        # Nj=1.4687, Qy=-0.71755, Qz=-0.48234, T=-36.432, My=-84.753, Mz=96.294

        model.statics_reactions(m)

        # for jid in [3, 4]:
        #     j = m["joints"][jid]
        #     print(f"Joint {jid}:")
        #     print(
        #         f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, Rz={j['reactions'][2]:.5}, Mx={j['reactions'][3]:.5}, My={j['reactions'][4]:.5}, Mz={j['reactions'][5]:.5}: "
        #     )
        j = m["joints"][3]
        self.assertAlmostEqual(j["reactions"][0], -1.1039, places=3)
        self.assertAlmostEqual(j["reactions"][1], -0.43381, places=3)
        self.assertAlmostEqual(j["reactions"][2], 0.21393, places=3)
        self.assertAlmostEqual(j["reactions"][3], 48.333, places=3)
        self.assertAlmostEqual(j["reactions"][4], -17.656, places=3)
        self.assertAlmostEqual(j["reactions"][5], 96.284, places=3)

        # plots.plot_setup(m)
        # plots.plot_members(m)
        # plots.plot_beam_orientation(m, 20)
        # plots.plot_deformations(m, 80.0)
        # # ax = plots.plot_shear_forces(m, scale=0.50e-3)
        # # ax.set_title('Shear forces')
        # plots.show(m)

    def test_linked_cantilevers_free(self):
        # SI units
        L = 3.0
        H = 0.3
        B = 0.2
        E = 200e9
        G = 80e9
        P = 60000

        A, Ix, Iy, Iz, J = section.rectangle(H, B)
        xz_vector = [0, 0, 1]
        sect_1 = section.beam_3d_section(
            "sect_1", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )
        A, Ix, Iy, Iz, J = section.rectangle(H, B / 2)
        xz_vector = [0, 0, 1]
        sect_2 = section.beam_3d_section(
            "sect_2", E=E, G=G, A=A, Ix=Ix, Iy=Iy, Iz=Iz, J=J, xz_vector=xz_vector
        )

        m = model.create(3)

        model.add_joint(m, 1, [0.0, 0.0, 0.0])
        model.add_joint(m, 2, [L, 0, 0.0])
        model.add_joint(m, 3, [2 * L, 0, 0])
        model.add_joint(m, 4, [L, 0, 0.0])

        model.add_support(m["joints"][1], model.CLAMPED)
        model.add_support(m["joints"][3], model.CLAMPED)

        model.add_beam_member(m, 1, [1, 2], sect_1)
        model.add_beam_member(m, 2, [3, 4], sect_2)

        model.add_load(m["joints"][4], model.U3, -P)

        model.add_link(m["joints"][2], m["joints"][4], model.U1)
        model.add_link(m["joints"][2], m["joints"][4], model.U2)
        model.add_link(m["joints"][2], m["joints"][4], model.U3)

        model.number_dofs(m)

        print("Number of free degrees of freedom = ", m["nfreedof"])
        print("Number of all degrees of freedom = ", m["ntotaldof"])

        # print([j['dof'] for j in m['joints'].values()])

        model.solve(m)

        for id in [2, 4]:
            j = m["joints"][id]
            print(id, j["displacements"])

        for k in m["beam_members"].keys():
            member = m["beam_members"][k]
            connectivity = member["connectivity"]
            i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
            f = beam.beam_3d_end_forces(member, i, j)
            print(f"Member {k}: ")
            print(
                f" Joint {connectivity[0]}: N={f['Ni']:.5}, Qy={f['Qyi']:.5}, Qz={f['Qzi']:.5}, T={f['Ti']:.5}, My={f['Myi']:.5}, Mz={f['Mzi']:.5}: "
            )
            print(
                f" Joint {connectivity[1]}: N={f['Nj']:.5}, Qy={f['Qyj']:.5}, Qz={f['Qzj']:.5}, T={f['Tj']:.5}, My={f['Myj']:.5}, Mz={f['Mzj']:.5}: "
            )

        model.statics_reactions(m)

        for jid in [1, 3]:
            j = m["joints"][jid]
            print(f"Joint {jid}:")
            print(
                f"   Rx={j['reactions'][0]:.5}, Ry={j['reactions'][1]:.5}, Rz={j['reactions'][2]:.5}, Mx={j['reactions'][3]:.5}, My={j['reactions'][4]:.5}, Mz={j['reactions'][5]:.5}: "
            )
        j1 = m["joints"][1]
        j3 = m["joints"][3]
        self.assertAlmostEqual(j1["reactions"][2] + j3["reactions"][2], P)

        # plots.plot_setup(m)
        # plots.plot_members(m)
        # plots.plot_member_numbers(m)
        # plots.plot_deformations(m, 100.0)
        # plots.plot_beam_orientation(m, 0.5)
        # # plots.plot_moments(m, 0.00001, "y")
        # # plots.plot_moments(m, 0.00001, "z")
        # # ax = plots.plot_shear_forces(m, scale=0.50e-3)
        # # ax.set_title('Shear forces')
        # plots.show(m)


def main():
    unittest.main()


if __name__ == "__main__":
    main()
