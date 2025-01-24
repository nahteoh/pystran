"""
Implement simple plots for truss and beam structures.
"""

import matplotlib.pyplot as plt
from numpy import linspace, dot, zeros
from pystran.beam import (
    beam_2d_member_geometry,
    beam_2d_shape_fun,
    beam_2d_moment,
    beam_2d_shear_force,
)
from pystran.beam import (
    beam_3d_member_geometry,
    beam_3d_xz_shape_fun,
    beam_3d_xy_shape_fun,
    beam_3d_moment,
    beam_3d_shear_force,
)


def plot_setup(m):
    """
    Setup the plot.
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
        ax.plot(xs, ys, "b-")
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
        ax.plot(xs, ys, zs, "b-")
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
        V = beam_2d_shear_force(member, i, j, xi)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs = zeros(2)
        ys = zeros(2)
        xs[0] = x[0]
        xs[1] = x[0] + scale * V * e_z[0]
        ys[0] = x[1]
        ys[1] = x[1] + scale * V * e_z[1]
        ax.plot(xs, ys, "b-")
        if xi == 0.0:
            ax.text(xs[1], ys[1], str(f"{V[0]:.5}"))
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
        ax.plot(xs, ys, zs, "b-")
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


def show(m):
    """
    Show the plot.
    """
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    if m["dim"] == 3:
        ax.set_zlabel("Z")
    plt.show()
