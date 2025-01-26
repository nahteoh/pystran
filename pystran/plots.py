"""
Implement simple plots for truss and beam structures.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Arc, RegularPolygon
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy
from numpy import linspace, dot, zeros, array
from numpy import radians as rad
from numpy.linalg import norm
from pystran.model import U1, U2, U3, UR1, UR2, UR3
from pystran.model import ALL_DOFS
from pystran.model import TRANSLATION_DOFS
from pystran.model import ndof_per_joint
from pystran.truss import (
    truss_member_geometry,
    truss_strain_displacement,
    truss_axial_force,
)
from pystran.beam import (
    beam_2d_member_geometry,
    beam_2d_shape_fun,
    beam_2d_moment,
    beam_2d_shear_force,
    beam_2d_axial_force,
)
from pystran.beam import (
    beam_3d_member_geometry,
    beam_3d_xz_shape_fun,
    beam_3d_xy_shape_fun,
    beam_3d_moment,
    beam_3d_shear_force,
    beam_3d_torsion_moment,
    beam_3d_axial_force,
)


# fig = plt.figure(figsize=(9,9))
# ax = plt.gca()
def _drawCirc(ax, radius, centX, centY, angle_, theta2_, sense, color_="black"):
    # ========Line
    arc = Arc(
        [centX, centY],
        radius,
        radius,
        angle=angle_,
        theta1=0,
        theta2=theta2_,
        capstyle="round",
        linestyle="-",
        lw=2,
        color=color_,
    )
    ax.add_patch(arc)

    # ========Create the arrow head
    if sense > 0:
        endX = centX + (radius / 2) * numpy.cos(
            rad(theta2_ + angle_)
        )  # Do trig to determine end position
        endY = centY + (radius / 2) * numpy.sin(rad(theta2_ + angle_))
    else:
        endX = centX + (radius / 2) * numpy.cos(
            rad(angle_)
        )  # Do trig to determine end position
        endY = centY + (radius / 2) * numpy.sin(rad(angle_))

    ax.add_patch(  # Create triangle as arrow head
        RegularPolygon(
            (endX, endY),  # (x,y)
            3,  # number of vertices
            radius=radius / 9,  # radius
            orientation=rad(angle_ + theta2_),  # orientation
            color=color_,
        )
    )
    # ax.set_xlim([centX - radius, centY + radius]) and ax.set_ylim(
    #     [centY - radius, centY + radius]
    # )
    # Make sure you keep the axes scaled or else arrow will distort


# _drawCirc(ax,1,1,1,0,250)
# _drawCirc(ax,2,1,1,90,330,color_='blue')
# plt.show()


# From: https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c
class _Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)
        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)
        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return numpy.min(zs)


def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    """Add an 3d arrow to an `Axes3D` instance."""

    arrow = _Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


setattr(Axes3D, "arrow3D", _arrow3D)


def plot_setup(m):
    """
    Setup the plot.

    This function creates a figure and an axis object. The axes are returned.
    """
    fig = plt.figure()
    if m["dim"] == 3:
        ax = fig.add_subplot(projection="3d")
    else:
        ax = fig.gca()
    return ax


def plot_members(m):
    """
    Plot the members of the structure.
    """
    ax = plt.gca()
    if m["dim"] == 3:
        for member in m["truss_members"].values():
            connectivity = member["connectivity"]
            i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
            ci, cj = i["coordinates"], j["coordinates"]
            plt.plot([ci[0], cj[0]], [ci[1], cj[1]], [ci[2], cj[2]], "k-")
        for member in m["beam_members"].values():
            connectivity = member["connectivity"]
            i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
            ci, cj = i["coordinates"], j["coordinates"]
            plt.plot([ci[0], cj[0]], [ci[1], cj[1]], [ci[2], cj[2]], "k-")
    else:
        for member in m["truss_members"].values():
            connectivity = member["connectivity"]
            i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
            ci, cj = i["coordinates"], j["coordinates"]
            plt.plot([ci[0], cj[0]], [ci[1], cj[1]], "k-")
        for member in m["beam_members"].values():
            connectivity = member["connectivity"]
            i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
            ci, cj = i["coordinates"], j["coordinates"]
            plt.plot([ci[0], cj[0]], [ci[1], cj[1]], "k-")
    return ax


def _plot_2d_beam_deflection(ax, member, i, j, scale):
    di, dj = i["displacements"], j["displacements"]
    ci, cj = i["coordinates"], j["coordinates"]
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    wi = dot(di[0:2], e_z)
    thi = di[2]
    wj = dot(dj[0:2], e_z)
    thj = dj[2]
    n = 20
    xs = zeros(n)
    ys = zeros(n)
    for s, xi in enumerate(linspace(-1, +1, n)):
        N = beam_2d_shape_fun(xi)
        w = N[0] * wi + (h / 2) * N[1] * thi + N[2] * wj + (h / 2) * N[3] * thj
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs[s] = x[0] + scale * w * e_z[0]
        ys[s] = x[1] + scale * w * e_z[1]
    ax.plot(xs, ys, "m-")


def _plot_3d_beam_deflection(ax, member, i, j, scale):
    sect = member["section"]
    di, dj = i["displacements"], j["displacements"]
    ci, cj = i["coordinates"], j["coordinates"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    ui = dot(di[0:3], e_x)
    uj = dot(dj[0:3], e_x)
    wi = dot(di[0:3], e_z)
    thyi = dot(di[3:6], e_y)
    wj = dot(dj[0:3], e_z)
    thyj = dot(dj[3:6], e_y)
    vi = dot(di[0:3], e_y)
    thzi = dot(di[3:6], e_z)
    vj = dot(dj[0:3], e_y)
    thzj = dot(dj[3:6], e_z)
    n = 20
    xs = zeros(n)
    ys = zeros(n)
    zs = zeros(n)
    for s, xi in enumerate(linspace(-1, +1, n)):
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs[s] = x[0]
        ys[s] = x[1]
        zs[s] = x[2]
        xs[s] += scale * ((1 - xi) / 2 * ui + (1 + xi) / 2 * uj) * e_x[0]
        ys[s] += scale * ((1 - xi) / 2 * ui + (1 + xi) / 2 * uj) * e_x[1]
        zs[s] += scale * ((1 - xi) / 2 * ui + (1 + xi) / 2 * uj) * e_x[2]
        N = beam_3d_xz_shape_fun(xi)
        w = N[0] * wi + (h / 2) * N[1] * thyi + N[2] * wj + (h / 2) * N[3] * thyj
        xs[s] += scale * w * e_z[0]
        ys[s] += scale * w * e_z[1]
        zs[s] += scale * w * e_z[2]
        N = beam_3d_xy_shape_fun(xi)
        v = N[0] * vi + (h / 2) * N[1] * thzi + N[2] * vj + (h / 2) * N[3] * thzj
        xs[s] += scale * v * e_y[0]
        ys[s] += scale * v * e_y[1]
        zs[s] += scale * v * e_y[2]
    ax.plot(xs, ys, zs, "m-")


def plot_deformations(m, scale=1.0):
    """
    Plot the deformations in the truss and beam members.
    """
    ax = plt.gca()
    for member in m["truss_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        di, dj = i["displacements"], j["displacements"]
        ci, cj = i["coordinates"], j["coordinates"]
        if m["dim"] == 3:
            ax.plot(
                [ci[0] + scale * di[0], cj[0] + scale * dj[0]],
                [ci[1] + scale * di[1], cj[1] + scale * dj[1]],
                [ci[2] + scale * di[2], cj[2] + scale * dj[2]],
                "m-",
            )
        else:
            ax.plot(
                [ci[0] + scale * di[0], cj[0] + scale * dj[0]],
                [ci[1] + scale * di[1], cj[1] + scale * dj[1]],
                "m-",
            )
    for member in m["beam_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        if m["dim"] == 3:
            _plot_3d_beam_deflection(ax, member, i, j, scale)
        else:
            _plot_2d_beam_deflection(ax, member, i, j, scale)
    return ax


def _plot_member_numbers_2d(m):
    ax = plt.gca()
    for jid in m["truss_members"].keys():
        member = m["truss_members"][jid]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        ci, cj = i["coordinates"], j["coordinates"]
        xm = (ci + cj) / 2.0
        ax.text(xm[0], xm[1], str(jid))
    for jid in m["beam_members"].keys():
        member = m["beam_members"][jid]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        ci, cj = i["coordinates"], j["coordinates"]
        xm = (ci + cj) / 2.0
        ax.text(xm[0], xm[1], str(jid))
    return ax


def _plot_member_numbers_3d(m):
    ax = plt.gca()
    for jid in m["truss_members"].keys():
        member = m["truss_members"][jid]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        ci, cj = i["coordinates"], j["coordinates"]
        xm = (ci + cj) / 2.0
        ax.text(xm[0], xm[1], xm[2], str(jid), "z")
    for jid in m["beam_members"].keys():
        member = m["beam_members"][jid]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        ci, cj = i["coordinates"], j["coordinates"]
        xm = (ci + cj) / 2.0
        ax.text(xm[0], xm[1], xm[2], str(jid), "z")
    return ax


def plot_member_numbers(m):
    """
    Plot the member numbers.
    """
    if m["dim"] == 3:
        ax = _plot_member_numbers_3d(m)
    else:
        ax = _plot_member_numbers_2d(m)
    return ax


def plot_joint_numbers(m):
    """
    Plot the joint numbers.
    """
    ax = plt.gca()
    for j in m["joints"].values():
        if m["dim"] == 3:
            ax.plot(j["coordinates"][0], j["coordinates"][1], j["coordinates"][2], "ro")
            ax.text(
                j["coordinates"][0],
                j["coordinates"][1],
                j["coordinates"][2],
                str(j["jid"]),
            )
        else:
            ax.plot(j["coordinates"][0], j["coordinates"][1], "ro")
            ax.text(j["coordinates"][0], j["coordinates"][1], str(j["jid"]))
    return ax


def _plot_2d_beam_moments(ax, member, i, j, scale):
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    ci, cj = i["coordinates"], j["coordinates"]
    n = 13
    for s, xi in enumerate(linspace(-1, +1, n)):
        M = beam_2d_moment(member, i, j, xi)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        # The convention: moment is plotted next to fibers in tension
        xs = zeros(2)
        ys = zeros(2)
        xs[0] = x[0]
        xs[1] = x[0] + scale * M * e_z[0]
        ys[0] = x[1]
        ys[1] = x[1] + scale * M * e_z[1]
        ax.plot(xs, ys, "r-" if (M > 0) else "b-")
        if xi == -1.0:
            ax.text(xs[1], ys[1], str(f"{M[0]:.5}"))
        elif xi == +1.0:
            ax.text(xs[1], ys[1], str(f"{M[0]:.5}"))
    return ax


def _plot_3d_beam_moments(ax, member, i, j, axis, scale):
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    ci, cj = i["coordinates"], j["coordinates"]
    n = 13
    dirv = e_y
    if axis == "y":
        dirv = e_z
    for s, xi in enumerate(linspace(-1, +1, n)):
        M = beam_3d_moment(member, i, j, axis, xi)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        zs = zeros(2)
        xs[0] = x[0]
        xs[1] = xs[0] + scale * M * dirv[0]
        ys[0] = x[1]
        ys[1] = ys[0] + scale * M * dirv[1]
        zs[0] = x[2]
        zs[1] = zs[0] + scale * M * dirv[2]
        ax.plot(xs, ys, zs, "r-" if (M > 0) else "b-")
        if xi == -1.0:
            ax.text(xs[1], ys[1], zs[1], str(f"{M[0]:.5}"))
        elif xi == +1.0:
            ax.text(xs[1], ys[1], zs[1], str(f"{M[0]:.5}"))
    return ax


def plot_moments(m, scale=1.0, axis="y"):
    """
    Plot the moments in the beam members.

    Optional: axis = "y" or "z" (default is "y", which is suitable for 2d beams).
    """
    ax = plt.gca()
    for member in m["beam_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        if m["dim"] == 3:
            _plot_3d_beam_moments(ax, member, i, j, axis, scale)
        else:
            _plot_2d_beam_moments(ax, member, i, j, scale)
    return ax


def _plot_2d_beam_shear_forces(ax, member, i, j, scale):
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    ci, cj = i["coordinates"], j["coordinates"]
    n = 13
    for s, xi in enumerate(linspace(-1, +1, n)):
        Q = beam_2d_shear_force(member, i, j, xi)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        xs[0] = x[0]
        xs[1] = x[0] + scale * Q * e_z[0]
        ys[0] = x[1]
        ys[1] = x[1] + scale * Q * e_z[1]
        ax.plot(xs, ys, "r-" if (Q > 0) else "b-")
        if xi == 0.0:
            ax.text(xs[1], ys[1], str(f"{Q[0]:.5}"))
    return ax


def _plot_3d_beam_shear_forces(ax, member, i, j, axis, scale):
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    ci, cj = i["coordinates"], j["coordinates"]
    n = 13
    dirv = e_z
    if axis == "y":
        dirv = e_y
    for s, xi in enumerate(linspace(-1, +1, n)):
        Q = beam_3d_shear_force(member, i, j, axis, xi)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        zs = zeros(2)
        xs[0] = x[0]
        xs[1] = xs[0] + scale * Q * dirv[0]
        ys[0] = x[1]
        ys[1] = ys[0] + scale * Q * dirv[1]
        zs[0] = x[2]
        zs[1] = zs[0] + scale * Q * dirv[2]
        ax.plot(xs, ys, zs, "r-" if (Q > 0) else "b-")
        if xi == -1.0:
            ax.text(xs[1], ys[1], zs[1], str(f"{Q[0]:.5}"))
        elif xi == +1.0:
            ax.text(xs[1], ys[1], zs[1], str(f"{Q[0]:.5}"))
    return ax


def plot_shear_forces(m, scale=1.0, axis="z"):
    """
    Plot the shear forces in the beam members.

    Optional: axis = "y" or "z" (default is "z", which is suitable for 2d beams).
    """
    ax = plt.gca()
    for member in m["beam_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        if m["dim"] == 3:
            _plot_3d_beam_shear_forces(ax, member, i, j, axis, scale)
        else:
            _plot_2d_beam_shear_forces(ax, member, i, j, scale)
    return ax


def _plot_2d_beam_axial_forces(ax, member, i, j, scale):
    e_x, e_z, h = beam_2d_member_geometry(i, j)
    ci, cj = i["coordinates"], j["coordinates"]
    n = 13
    for s, xi in enumerate(linspace(-1, +1, n)):
        N = beam_2d_axial_force(member, i, j)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        xs[0] = x[0]
        xs[1] = x[0] + scale * N * e_z[0]
        ys[0] = x[1]
        ys[1] = x[1] + scale * N * e_z[1]
        ax.plot(xs, ys, "r-" if (N > 0) else "b-")
        if xi == 0.0:
            ax.text(xs[1], ys[1], str(f"{N[0]:.5}"))
    return ax


def _plot_2d_truss_axial_forces(ax, member, i, j, scale):
    e_x, h = truss_member_geometry(i, j)
    ci, cj = i["coordinates"], j["coordinates"]
    N = truss_axial_force(member, i, j)
    e_z = array([e_x[1], -e_x[0]])
    n = 13
    for s, xi in enumerate(linspace(-1, +1, n)):
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        xs[0] = x[0]
        xs[1] = x[0] + scale * N * e_z[0]
        ys[0] = x[1]
        ys[1] = x[1] + scale * N * e_z[1]
        ax.plot(xs, ys, "r-" if (N > 0) else "b-")
        if xi == 0.0:
            ax.text(xs[1], ys[1], str(f"{N[0]:.5}"))
    return ax


def _plot_3d_beam_axial_forces(ax, member, i, j, scale):
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    ci, cj = i["coordinates"], j["coordinates"]
    n = 13
    dirv = e_z
    for s, xi in enumerate(linspace(-1, +1, n)):
        N = beam_3d_axial_force(member, i, j)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        zs = zeros(2)
        xs[0] = x[0]
        xs[1] = xs[0] + scale * N * dirv[0]
        ys[0] = x[1]
        ys[1] = ys[0] + scale * N * dirv[1]
        zs[0] = x[2]
        zs[1] = zs[0] + scale * N * dirv[2]
        ax.plot(xs, ys, zs, "r-" if (N > 0) else "b-")
        if xi == 0.0:
            ax.text(xs[1], ys[1], zs[1], str(f"{N[0]:.5}"))
    return ax


def plot_axial_forces(m, scale=1.0):
    """
    Plot the axial forces in the members.

    Optional: axis = "y" or "z" (default is "z", which is suitable for 2d beams).
    """
    ax = plt.gca()
    for member in m["truss_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        if m["dim"] == 3:
            _plot_3d_truss_axial_forces(ax, member, i, j, scale)
        else:
            _plot_2d_truss_axial_forces(ax, member, i, j, scale)
    for member in m["beam_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        if m["dim"] == 3:
            _plot_3d_beam_axial_forces(ax, member, i, j, scale)
        else:
            _plot_2d_beam_axial_forces(ax, member, i, j, scale)
    return ax


def _plot_3d_beam_torsion_moments(ax, member, i, j, scale):
    sect = member["section"]
    e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
    ci, cj = i["coordinates"], j["coordinates"]
    n = 13
    dirv = e_z
    for s, xi in enumerate(linspace(-1, +1, n)):
        T = beam_3d_torsion_moment(member, i, j)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        zs = zeros(2)
        xs[0] = x[0]
        xs[1] = xs[0] + scale * T * dirv[0]
        ys[0] = x[1]
        ys[1] = ys[0] + scale * T * dirv[1]
        zs[0] = x[2]
        zs[1] = zs[0] + scale * T * dirv[2]
        ax.plot(xs, ys, zs, "r-" if (T > 0) else "b-")
        if xi == 0.0:
            ax.text(xs[1], ys[1], zs[1], str(f"{T[0]:.5}"))
    return ax


def plot_torsion_moments(m, scale=1.0):
    """
    Plot the torsion moments in the 3D beam members.
    """
    ax = plt.gca()
    for member in m["beam_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        if m["dim"] == 3:
            _plot_3d_beam_torsion_moments(ax, member, i, j, scale)
    return ax


def plot_beam_orientation(m, scale=1.0):
    """
    Plot the beam orientations as cartesian triplets.
    """
    ax = plt.gca()
    for member in m["beam_members"].values():
        sect = member["section"]
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        ci, cj = i["coordinates"], j["coordinates"]
        xm = (ci + cj) / 2.0
        if m["dim"] == 3:
            e_x, e_y, e_z, h = beam_3d_member_geometry(i, j, sect["xz_vector"])
            xs = zeros(2)
            ys = zeros(2)
            zs = zeros(2)
            xs[0] = xm[0]
            ys[0] = xm[1]
            zs[0] = xm[2]
            xs[1] = xs[0] + scale * e_x[0]
            ys[1] = ys[0] + scale * e_x[1]
            zs[1] = zs[0] + scale * e_x[2]
            ax.plot(xs, ys, zs, "r-", lw=3)
            xs = zeros(2)
            ys = zeros(2)
            zs = zeros(2)
            xs[0] = xm[0]
            ys[0] = xm[1]
            zs[0] = xm[2]
            xs[1] = xs[0] + scale * e_y[0]
            ys[1] = ys[0] + scale * e_y[1]
            zs[1] = zs[0] + scale * e_y[2]
            ax.plot(xs, ys, zs, "g-", lw=3)
            xs = zeros(2)
            ys = zeros(2)
            zs = zeros(2)
            xs[0] = xm[0]
            ys[0] = xm[1]
            zs[0] = xm[2]
            xs[1] = xs[0] + scale * e_z[0]
            ys[1] = ys[0] + scale * e_z[1]
            zs[1] = zs[0] + scale * e_z[2]
            ax.plot(xs, ys, zs, "b-", lw=3)
        else:
            e_x, e_z, h = beam_2d_member_geometry(i, j)
            xs = zeros(2)
            ys = zeros(2)
            xs[0] = xm[0]
            ys[0] = xm[1]
            xs[1] = xs[0] + scale * e_x[0]
            ys[1] = ys[0] + scale * e_x[1]
            ax.plot(xs, ys, "r-", lw=3)
            xs = zeros(2)
            ys = zeros(2)
            xs[0] = xm[0]
            ys[0] = xm[1]
            xs[1] = xs[0] + scale * e_z[0]
            ys[1] = ys[0] + scale * e_z[1]
            ax.plot(xs, ys, "b-", lw=3)
    return ax


def plot_applied_forces(m, scale=1.0):
    """
    Plot the applied forces at the joints.

    - `m` = model dictionary,
    - `scale` = scale factor for the arrows.
    """
    ax = plt.gca()
    dim = m["dim"]
    ndpn = ndof_per_joint(m)
    for j in m["joints"].values():
        if "loads" in j and j["loads"]:
            F = zeros((dim,))
            for d in j["loads"].keys():
                if d < dim:
                    F[d] = j["loads"][d]
                if norm(F) > 0:
                    if dim == 2:
                        x, y = j["coordinates"]
                        u, v = F
                        al = scale * norm(F)
                        ax.arrow(
                            x,
                            y,
                            scale * u,
                            scale * v,
                            head_width=al / 5,
                            head_length=al / 5,
                            color="cyan",
                        )
                    else:
                        x, y, z = j["coordinates"]
                        u, v, w = F
                        ax.arrow3D(
                            x,
                            y,
                            z,
                            scale * u,
                            scale * v,
                            scale * w,
                            mutation_scale=20,
                            arrowstyle="-|>",
                            color="cyan",
                        )
    return ax


def plot_applied_moments(m, scale=1.0, radius=0.1):
    """
    Plot the applied moments at the joints.

    - `m` = model dictionary,
    - `scale` = scale factor for the arrows. Moments are rendered with double
      arrows.

    Optional: `radius` = radius of the circle (2D only).
    """
    ax = plt.gca()
    dim = m["dim"]
    ndpn = ndof_per_joint(m)
    for j in m["joints"].values():
        if "loads" in j and j["loads"]:
            M = zeros((ndpn - dim,))
            for d in j["loads"].keys():
                if d >= dim:
                    M[d - dim] = j["loads"][d]
                if norm(M) > 0:
                    if dim == 2:
                        x, y = j["coordinates"]
                        if M > 0:
                            st = -110
                            dl = 210
                            sense = +1
                        else:
                            st = 80
                            dl = 210
                            sense = -1
                        _drawCirc(ax, radius, x, y, st, dl, sense, color_="cyan")
                    else:
                        x, y, z = j["coordinates"]
                        u, v, w = M
                        ax.arrow3D(
                            x,
                            y,
                            z,
                            scale * u,
                            scale * v,
                            scale * w,
                            mutation_scale=20,
                            arrowstyle="-|>",
                            color="cyan",
                        )
                        ax.arrow3D(
                            x,
                            y,
                            z,
                            scale * 0.9 * u,
                            scale * 0.9 * v,
                            scale * 0.9 * w,
                            mutation_scale=20,
                            arrowstyle="-|>",
                            color="cyan",
                        )
    return ax


def plot_translation_supports(m, scale=1.0, shortest_arrow=1.0e-6):
    """
    Plot the translation supports at the joints.

    - `m` = model dictionary,
    - `scale` = scale factor for the arrows.
    """
    ax = plt.gca()
    dim = m["dim"]
    ndpn = ndof_per_joint(m)
    for j in m["joints"].values():
        if "supports" in j and j["supports"]:
            for d in j["supports"].keys():
                U = zeros((dim,))
                if d < dim:
                    v = j["supports"][d]
                    if v == 0:
                        v = shortest_arrow
                    U[d] = v
                if dim == 2:
                    x, y = j["coordinates"]
                    u, v = F
                    al = scale * norm(F)
                    ax.arrow(
                        x,
                        y,
                        scale * u,
                        scale * v,
                        head_width=al / 5,
                        head_length=al / 5,
                        color="cyan",
                    )
                else:
                    x, y, z = j["coordinates"]
                    u, v, w = U
                    ax.arrow3D(
                        x,
                        y,
                        z,
                        scale * u,
                        scale * v,
                        scale * w,
                        mutation_scale=20,
                        arrowstyle="-|>",
                        color="cyan",
                    )
    return ax


def plot_rotation_supports(m, scale=1.0, radius=0.1, shortest_arrow=1.0e-6):
    """
    Plot the rotation supports at the joints.

    - `m` = model dictionary,
    - `scale` = scale factor for the arrows. Moments are rendered with double
      arrows.

    Optional: `radius` = radius of the circle (2D only).
    """
    ax = plt.gca()
    dim = m["dim"]
    ndpn = ndof_per_joint(m)
    if ndpn == dim:
        return
    for j in m["joints"].values():
        if "supports" in j and j["supports"]:
            for d in j["supports"].keys():
                R = zeros((ndpn - dim,))
                if d >= dim:
                    v = j["supports"][d]
                    if v == 0:
                        v = shortest_arrow
                    R[d - dim] = v
                    if dim == 2:
                        x, y = j["coordinates"]
                        if R >= 0:
                            st = -110
                            dl = 210
                            sense = +1
                        else:
                            st = 80
                            dl = 210
                            sense = -1
                        _drawCirc(ax, radius, x, y, st, dl, sense, color_="cyan")
                    else:
                        x, y, z = j["coordinates"]
                        u, v, w = scale * R
                        ax.arrow3D(
                            x,
                            y,
                            z,
                            u,
                            v,
                            w,
                            mutation_scale=20,
                            arrowstyle="-|>",
                            color="cyan",
                        )
                        if norm(R) > 0:
                            ax.arrow3D(
                                x,
                                y,
                                z,
                                0.9 * u,
                                0.9 * v,
                                0.9 * w,
                                mutation_scale=20,
                                arrowstyle="-|>",
                                color="cyan",
                            )
    return ax


def show(m):
    """
    Show the plot.
    """
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    if m["dim"] == 3:
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
    plt.show()
