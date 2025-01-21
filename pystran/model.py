"""
Define the functions for defining and manipulating a model.
"""

from math import sqrt
import numpy
from numpy import array, zeros, dot
import scipy
import pystran.section

# These are the designations of degrees of freedom (directions in space), plus a
# special designation for a clamped support.
U1 = 0
U2 = 1
U3 = 2
UR1 = 3
UR2 = 4
UR3 = 5
CLAMPED = 100
PINNED = 200


def create(dim=2):
    """
    Create a new model.

    Supply the dimension of the model (2 or 3).
    """
    m = dict()
    m["dim"] = dim  # Dimension of the model
    m["joints"] = dict()
    m["truss_members"] = dict()
    m["beam_members"] = dict()

    global U1
    global U2
    global U3
    global UR1
    global UR2
    global UR3
    if m["dim"] == 2:
        U1 = 0
        U2 = 1
        UR3 = 2
        U3 = -1000  # invalid
        UR1 = -1000  # invalid
        UR2 = -1000  # invalid
    else:
        U1 = 0
        U2 = 1
        U3 = 2
        UR1 = 3
        UR2 = 4
        UR3 = 5
    return m


def add_joint(m, identifier, coordinates):
    """
    Add a joint to the model.
    """
    if identifier in m["joints"]:
        raise RuntimeError("Joint already exists")
    else:
        m["joints"][identifier] = {"coordinates": array(coordinates)}
    if m["joints"][identifier]["coordinates"].shape != (m["dim"],):
        raise RuntimeError("Coordinate dimension mismatch")
    return None


def add_truss_member(m, identifier, connectivity, sect):
    """
    Add a truss member to the model.
    """
    if identifier in m["truss_members"]:
        raise RuntimeError("Truss member already exists")
    else:
        m["truss_members"][identifier] = {
            "connectivity": array(connectivity, dtype=numpy.int32),
            "section": sect,
        }
    return None


def add_beam_member(m, identifier, connectivity, sect):
    """
    Add a beam member to the model.
    """
    if identifier in m["beam_members"]:
        raise RuntimeError("Beam member already exists")
    else:
        m["beam_members"][identifier] = {
            "connectivity": array(connectivity, dtype=numpy.int32),
            "section": sect,
        }
    return None


def add_support(j, dof, value=0.0):
    """
    Add a support to a joint.
    """
    if "supports" not in j:
        j["supports"] = dict()
    dim = len(j["coordinates"])
    if dof == CLAMPED:
        if dim == 2:
            j["supports"] = {U1: 0.0, U2: 0.0, UR3: 0.0}
        else:
            j["supports"] = {U1: 0.0, U2: 0.0, U3: 0.0, UR1: 0.0, UR2: 0.0, UR3: 0.0}
    elif dof == PINNED:
        if dim == 2:
            j["supports"] = {U1: 0.0, U2: 0.0}
        else:
            j["supports"] = {U1: 0.0, U2: 0.0, U3: 0.0}
    else:
        j["supports"][dof] = value
    return None


def add_load(j, dof, value):
    """
    Add a load to a joint.
    """
    if "loads" not in j:
        j["loads"] = dict()
    j["loads"][dof] = value
    return None


def add_mass(j, dof, value):
    """
    Add a mass to a joint.
    """
    if "mass" not in j:
        j["mass"] = dict()
    j["mass"][dof] = value
    return None


def number_dofs(m):
    """
    Number degrees of freedom.
    """
    # Determine the number of degrees of freedom per joint
    translation_only = not m["beam_members"]
    if translation_only:
        ndof_per_joint = m["dim"]
    else:
        if m["dim"] == 2:
            ndof_per_joint = 3
        else:
            ndof_per_joint = 6
    # Number the free degrees of freedom first
    n = 0
    for j in m["joints"].values():
        j["dof"] = array([i for i in range(ndof_per_joint)], dtype=numpy.int32)
        for d in range(ndof_per_joint):
            if "supports" not in j or d not in j["supports"]:
                j["dof"][d] = n
                n += 1
    m["nfreedof"] = n
    # Number all prescribed degrees of freedom
    for j in m["joints"].values():
        for d in range(ndof_per_joint):
            if "supports" in j and d in j["supports"]:
                j["dof"][d] = n
                n += 1
    m["ntotaldof"] = n
    return None


def solve(m):
    """
    The default solve procedure of the discrete model.
    """
    return solve_statics(m)


def solve_statics(m):
    """
    Solve the static equilibrium of the discrete model.
    """
    nt, nf = m["ntotaldof"], m["nfreedof"]
    # Assemble global stiffness matrix
    K = zeros((nt, nt))
    for member in m["truss_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        pystran.truss.assemble_stiffness(K, member, i, j)
    for member in m["beam_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        pystran.beam.assemble_stiffness(K, member, i, j)

    m["K"] = K

    # Apply boundary conditions
    F = zeros(m["ntotaldof"])
    for joint in m["joints"].values():
        if "loads" in joint:
            for dof, value in joint["loads"].items():
                gr = joint["dof"][dof]
                F[gr] += value

    m["F"] = F

    U = zeros(m["ntotaldof"])
    for joint in m["joints"].values():
        if "supports" in joint:
            for dof, value in joint["supports"].items():
                gr = joint["dof"][dof]
                U[gr] = value
    # # Solve for displacements
    U[0:nf] = numpy.linalg.solve(K[0:nf, 0:nf], F[0:nf] - dot(K[0:nf, nf:nt], U[nf:nt]))

    m["U"] = U

    # # Assign displacements back to joints
    for joint in m["joints"].values():
        joint["displacements"] = U[joint["dof"]]

    return None


def statics_reactions(m):
    """
    Compute the reactions in the static equilibrium of the discrete model.
    """
    nt, nf = m["ntotaldof"], m["nfreedof"]

    K = m["K"]
    F = m["F"]
    U = m["U"]

    # Compute reactions from the partitioned stiffness matrix and the
    # partitioned displacement vector
    # R = dot(K[nf:nt, 0:nf], U[0:nf]) + dot(K[nf:nt, nf:nt], U[nf:nt])
    R = dot(K, U)

    for joint in m["joints"].values():
        if "supports" in joint:
            reactions = dict()
            for dof, value in joint["supports"].items():
                gr = joint["dof"][dof]
                reactions[dof] = R[gr]
            joint["reactions"] = reactions

    return None


def solve_free_vibration(m):
    """
    Solve the free vibration of the discrete model.
    """
    nt, nf = m["ntotaldof"], m["nfreedof"]
    # Assemble global stiffness matrix and mass matrix
    K = zeros((nt, nt))
    M = zeros((nt, nt))
    for member in m["truss_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        pystran.truss.assemble_stiffness(K, member, i, j)
        pystran.truss.assemble_mass(M, member, i, j)
    for member in m["beam_members"].values():
        connectivity = member["connectivity"]
        i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
        pystran.beam.assemble_stiffness(K, member, i, j)
        pystran.beam.assemble_mass(M, member, i, j)
    for j in m["joints"].values():
        if "mass" in j:
            for dof, value in j["mass"].items():
                gr = j["dof"][dof]
                M[gr, gr] += value

    m["K"] = K
    m["M"] = M

    U = zeros(m["ntotaldof"])
    for joint in m["joints"].values():
        if "supports" in joint:
            for dof, value in joint["supports"].items():
                gr = joint["dof"][dof]
                U[gr] = 0.0
    m["U"] = U

    # Solved the eigenvalue problem
    eigvals, eigvecs = scipy.linalg.eigh(K[0:nf, 0:nf], M[0:nf, 0:nf])

    m["eigvals"] = eigvals
    m["frequencies"] = [sqrt(ev) / 2 / numpy.pi for ev in eigvals]
    m["eigvecs"] = eigvecs

    return


def copy_mode(m, mode):
    """
    Copy a mode to the displacement field of the model.
    """
    nf = m["nfreedof"]
    m["U"][0:nf] = m["eigvecs"][:, mode]
    for joint in m["joints"].values():
        joint["displacements"] = m["U"][joint["dof"]]
    return None
