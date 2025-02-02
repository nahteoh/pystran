"""
Define the functions for defining and manipulating a model.
"""

from math import sqrt
import numpy
from numpy import array, zeros, dot, mean, concatenate
import scipy
import pystran.section

U1 = 0
"""
This is a designation of the degree of freedom as translation along `X`.
"""
U2 = 1
"""
This is a designation of the degree of freedom as translation along `Z` (in 2D models)
or along `Y` (in 3D models).
"""
U3 = 2
"""
This is a designation of the degree of freedom as translation along `Z` (in 3D models).
"""
UR1 = 3
"""
This is a designation of the degree of freedom as rotation about  `X` (in 3D models).
"""
UR2 = 4
"""
This is a designation of the degree of freedom as rotation about  `Y` (in 3D models).
"""
UR3 = 5
"""
This is a designation of the degree of freedom as rotation about  `Y` (in 2D
models) or rotation about `Z`  (in 3D models).
"""
ALL_DOFS = 100
"""
This is a designation of all the degrees of freedom, translations and rotations
(`U1`, `U2`, `UR3` in 2D models, or `U1`, `U2`, `U3`, `UR1`, `UR2`, `UR3` in 3D
models). It may be used to specify the clamped condition for the joint.
"""
TRANSLATION_DOFS = 200
"""
This is a designation of the translation degrees of freedom (`U1`, `U2`,  in 2D
models, or `U1`, `U2`,  `U3` in 3D models). It may be used to specify the pinned
condition for the joint.
"""


def create(dim=2):
    """
    Create a new model.

    Supply the dimension of the model (2 or 3).
    """
    m = {}
    m["dim"] = dim  # Dimension of the model
    m["joints"] = {}
    m["truss_members"] = {}
    m["beam_members"] = {}

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


def add_joint(m, jid: int, coordinates):
    """
    Add a joint to the model.

    - `m` = the model,
    - `jid` = the joint identifier,
    - `coordinates` = the list of coordinates of the joint.
    """
    if jid in m["joints"]:
        raise RuntimeError("Joint already exists")
    coordinates = array(coordinates, dtype=numpy.float64)
    if coordinates.shape != (m["dim"],):
        raise RuntimeError("Coordinate dimension mismatch")
    m["joints"][jid] = {"jid": jid, "coordinates": array(coordinates)}


def add_truss_member(m, mid: int, connectivity, sect):
    """
    Add a truss member to the model.

    - `m` = the model,
    - `mid` = the member identifier (must be unique),
    - `connectivity` = the list of the joint identifiers,
    - `sect` = the section of the member.
    """
    if mid in m["truss_members"]:
        raise RuntimeError("Truss member already exists")
    m["truss_members"][mid] = {
        "connectivity": array(connectivity, dtype=numpy.int32),
        "section": sect,
    }


def add_beam_member(m, mid: int, connectivity, sect):
    """
    Add a beam member to the model.

    - `m` = the model,
    - `mid` = the member identifier (must be unique),
    - `connectivity` = the list of the joint identifiers,
    - `sect` = the section of the member.
    """
    if mid in m["beam_members"]:
        raise RuntimeError("Beam member already exists")
    else:
        m["beam_members"][mid] = {
            "connectivity": array(connectivity, dtype=numpy.int32),
            "section": sect,
        }


def _pinned_dofs(dim):
    if dim == 2:
        return [U1, U2]
    else:
        return [U1, U2, U3]


def _clamped_dofs(dim):
    if dim == 2:
        return [U1, U2, UR3]
    else:
        return [U1, U2, U3, UR1, UR2, UR3]


def _dofs_values(dim, dof, value):
    if dof == ALL_DOFS:
        return _clamped_dofs(dim), [0.0 for d in _clamped_dofs(dim)]
    elif dof == TRANSLATION_DOFS:
        return _pinned_dofs(dim), [0.0 for d in _pinned_dofs(dim)]
    return [dof], [value]


def add_support(j, dof, value=0.0):
    """
    Add a support to a joint.

    - `j` = the joint,
    - `dof` = the degree of freedom,
    - `value` = the amount of the support motion (default is zero).
    """
    if "supports" not in j:
        j["supports"] = {}
    dim = len(j["coordinates"])
    for d, v in zip(*_dofs_values(dim, dof, value)):
        j["supports"][d] = v


def add_load(j, dof, value):
    """
    Add a load to a joint.

    - `j` = the joint,
    - `dof` = the degree of freedom,
    - `value` = signed magnitude of the load.
    """
    if "loads" not in j:
        j["loads"] = {}
    j["loads"][dof] = value


def add_mass(j, dof, value):
    """
    Add a mass to a joint.

    - `j` = the joint,
    - `dof` = the degree of freedom,
    - `value` = magnitude of the mass.
    """
    if "mass" not in j:
        j["mass"] = {}
    j["mass"][dof] = value


def add_links(m, jids, dof):
    """
    Add links between all joints in the list `jids` in the direction `dof`.

    - `m` = the model,
    - `jids` = the list of joint identifiers,
    - `dof` = the degree of freedom at which the joints are to be linked.
    """
    # If any of the joints is supported, the link is not allowed at this stage
    # if "supports" in i or "supports" in j:
    #     raise RuntimeError("Link must be applied before supports")
    # Now add the mutual links between the joints
    for jid1 in jids:
        for jid2 in jids:
            if jid1 != jid2:
                j1 = m["joints"][jid1]
                if "links" not in j1:
                    j1["links"] = {}
                if jid2 not in j1["links"]:
                    j1["links"][jid2] = []
                for d, v in zip(*_dofs_values(m["dim"], dof, 0.0)):
                    j1["links"][jid2].append(d)


def bounding_box(m):
    """
    Compute the bounding box of the model.
    """
    dim = m["dim"]
    box = numpy.array(
        concatenate([[numpy.inf for i in range(dim)], [-numpy.inf for i in range(dim)]])
    )
    for j in m["joints"].values():
        cj = j["coordinates"]
        for i in range(len(cj)):
            box[i] = min(box[i], cj[i])
            box[i + dim] = max(box[i + dim], cj[i])
    return box


def characteristic_dimension(m):
    """
    Compute the characteristic dimension of the model.

    This is the average of the dimensions of the bounding box.
    """
    dim = m["dim"]
    box = bounding_box(m)
    dl = [box[i + dim] - box[i] for i in range(dim)]
    return mean(array(dl))


def _copy_dof_num_to_linked(m, j, d, n):
    if "links" in j:
        for k in j["links"].keys():
            o = m["joints"][k]
            if d in o["links"][j["jid"]]:
                o["dof"][d] = n


def _have_rotations(m):
    with_rotations = m["beam_members"]
    if with_rotations:
        return True
    for j in m["joints"].values():
        if "supports" in j and j["supports"]:
            for dof in j["supports"].keys():
                if dof == UR1 or dof == UR2 or dof == UR3:
                    return True
    return False


def ndof_per_joint(m):
    """
    How many degrees of freedom are there per joint?
    """
    ndpn = m["dim"]
    with_rotations = _have_rotations(m)
    if with_rotations:
        if m["dim"] == 2:
            ndpn = 3
        else:
            ndpn = 6
    return ndpn


def number_dofs(m):
    """
    Number degrees of freedom.
    """
    # Determine the number of degrees of freedom per joint
    ndpn = ndof_per_joint(m)
    # Generate arrays for storing the degrees of freedom
    for j in m["joints"].values():
        j["dof"] = zeros((ndpn,), dtype=numpy.int32) - 1
    # For each linked pair of joints, make sure they share the same supports
    for j in m["joints"].values():
        if "links" in j and "supports" in j:
            for k in j["links"].keys():
                o = m["joints"][k]
                if not "supports" in o:
                    o["supports"] = j["supports"].copy()
                if o["supports"] != j["supports"]:
                    raise RuntimeError("Linked joints must have the same supports")

    # Number the free degrees of freedom first
    n = 0
    for j in m["joints"].values():
        for d in range(ndpn):
            if ("supports" not in j) or (d not in j["supports"]):
                if j["dof"][d] < 0:
                    j["dof"][d] = n
                    _copy_dof_num_to_linked(m, j, d, n)
                    n += 1
    m["nfreedof"] = n
    # Number all prescribed degrees of freedom
    for j in m["joints"].values():
        for d in range(ndpn):
            if "supports" in j and d in j["supports"]:
                if j["dof"][d] < 0:
                    j["dof"][d] = n
                    _copy_dof_num_to_linked(m, j, d, n)
                    n += 1
    m["ntotaldof"] = n


def solve_statics(m):
    """
    Solve the static equilibrium of the discrete model.

    `number_dofs` must be called before this function.
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

    # Compute the active load vector
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
                if value != 0.0:
                    gr = joint["dof"][dof]
                    U[gr] = value
    # # Solve for displacements
    U[0:nf] = numpy.linalg.solve(K[0:nf, 0:nf], F[0:nf] - dot(K[0:nf, nf:nt], U[nf:nt]))

    m["U"] = U

    # # Assign displacements back to joints
    for joint in m["joints"].values():
        joint["displacements"] = U[joint["dof"]]


def statics_reactions(m):
    """
    Compute the reactions in the static equilibrium of the discrete model.

    The static solution must be obtained with `solve_statics` before calling
    this function.
    """
    K = m["K"]
    U = m["U"]

    # Compute reactions from the partitioned stiffness matrix and the
    # partitioned displacement vector
    # R = dot(K[nf:nt, 0:nf], U[0:nf]) + dot(K[nf:nt, nf:nt], U[nf:nt])
    R = dot(K, U)

    for joint in m["joints"].values():
        if "supports" in joint:
            reactions = {}
            for dof, value in joint["supports"].items():
                gr = joint["dof"][dof]
                reactions[dof] = R[gr]
            joint["reactions"] = reactions


def solve_free_vibration(m):
    """
    Solve the free vibration of the discrete model.

    The free vibration eigenvalue problem is solved for the eigenvalues and
    eigenvectors (can be retrieved as  `m["eigvals"]` and `m["eigvecs"]`). The
    frequencies are computed from the eigenvalues (can be retrieved as
    `m["frequencies"]`).

    `number_dofs` must be called before this function.
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


def set_solution(m, V):
    """
    Set the displacement solution from a vector. 

    - `m` = the model,
    - `V` = the displacement vector. Either of length for only the free degrees
      of freedom or for the total number of degrees of freedom.
    """
    nf = m["nfreedof"]
    nt = m["ntotaldof"]
    if len(V) == nf:
        m["U"][0:nf] = V
    elif len(V) == nt:
        m["U"][0:nf] = V
    else:
        raise RuntimeError("Invalid vector length")
    for joint in m["joints"].values():
        joint["displacements"] = m["U"][joint["dof"]]


def free_body_check(m):
    """
    Check the balance of the structure as a free body.

    All the active forces and moments together with the reactions at all the
    supports should sum to zero.

    `statics_reactions` must be called before this function as this calculation
    relies on the presence of reactions at the joints.
    """
    if m["dim"] == 2:
        nrbm = 3  # Number of rigid body modes: assume 2 translations, 1 rotation
        MZ = 2
        allforces = zeros(nrbm)
        for joint in m["joints"].values():
            c = joint["coordinates"]
            x, y = c[0], c[1]
            if "loads" in joint:
                for dof, value in joint["loads"].items():
                    if dof < MZ:
                        # Add contributions of forces to the moment
                        if dof == 0:
                            allforces[MZ] += -value * y
                        elif dof == 1:
                            allforces[MZ] += +value * x
                    else:
                        # Add contributions of forces and moments
                        allforces[dof] += value
            if "reactions" in joint:
                for dof, value in joint["reactions"].items():
                    if dof < MZ:
                        # Add contributions of forces to the moment
                        if dof == 0:
                            allforces[MZ] += -value * y
                        elif dof == 1:
                            allforces[MZ] += +value * x
                    else:
                        # Add contributions of forces and moments
                        allforces[dof] += value
        return allforces
    else:
        nrbm = 6  # Number of rigid body modes: assume 3 translations, 3 rotations
        MX, MY, MZ = 3, 4, 5
        allforces = zeros(nrbm)
        for joint in m["joints"].values():
            c = joint["coordinates"]
            x, y, z = c[0], c[1], c[2]
            if "loads" in joint:
                for dof, value in joint["loads"].items():
                    if dof < MX:
                        # Add contributions of forces to the moment
                        if dof == 0:
                            allforces[MY] += +value * z
                            allforces[MZ] += -value * y
                        elif dof == 1:
                            allforces[MX] += -value * z
                            allforces[MZ] += +value * x
                        else:
                            allforces[MY] += -value * x
                            allforces[MX] += +value * y
                    else:
                        # Add contributions of forces and moments
                        allforces[dof] += value
            if "reactions" in joint:
                for dof, value in joint["reactions"].items():
                    if dof < MX:
                        # Add contributions of forces to the moment
                        if dof == 0:
                            allforces[MY] += +value * z
                            allforces[MZ] += -value * y
                        elif dof == 1:
                            allforces[MX] += -value * z
                            allforces[MZ] += +value * x
                        else:
                            allforces[MY] += -value * x
                            allforces[MX] += +value * y
                    else:
                        # Add contributions of forces and moments
                        allforces[dof] += value
        return allforces


def refine_member(m, mid, n):
    """
    Refine a beam member by replacing it with `n` new members.

    The new joints are numbered starting from the maximum joint identifier in
    the model.

    - `m` = the model,
    - `mid` = the member identifier,
    - `n` = the number of new beam members to replace the old member with.
    """
    if n < 2:
        raise RuntimeError("Number of new members must be at least 2")
    member = m["beam_members"][mid]
    connectivity = member["connectivity"]
    i, j = m["joints"][connectivity[0]], m["joints"][connectivity[1]]
    ci, cj = i["coordinates"], j["coordinates"]
    # First replacement member
    start = -1.0 + 2.0 / n
    c = (-1 + start) / (-2) * ci + (1 + start) / 2 * cj
    newjid = max(m["joints"].keys()) + 1
    add_joint(m, newjid, c)
    newmid = max(m["beam_members"].keys()) + 1
    add_beam_member(m, newmid, [i["jid"], newjid], member["section"])
    prevjid = newjid
    for k in range(n - 2):
        start += 2.0 / n
        c = (-1 + start) / (-2) * ci + (1 + start) / 2 * cj
        newjid = max(m["joints"].keys()) + 1
        add_joint(m, newjid, c)
        newmid = max(m["beam_members"].keys()) + 1
        add_beam_member(m, newmid, [prevjid, newjid], member["section"])
        prevjid = newjid
    # Last replacement member
    newmid = max(m["beam_members"].keys()) + 1
    add_beam_member(m, newmid, [newjid, j["jid"]], member["section"])
    # Remove the old member
    del m["beam_members"][mid]
