"""
Define rigid link mechanical quantities.
"""

from numpy import reshape, outer, concatenate, zeros, dot, array, eye
from pystran import geometry
from pystran import assemble
from pystran import gauss


def rigid_link_stiffness(e_x, h, Gamma):
    r"""
    Compute rigid link stiffness matrix.

    The stiffness matrix is computed as

    $$
    K =\begin{bmatrix} 
          C^T\Gamma C & -C^T\Gamma  \\
          -\Gamma C & \Gamma \\
       \end{bmatrix}.
    $$

    Here $C$ is a matrix computed from the vector    $r = h e_x$, 
    which is the difference between the location of 
    the subordinate and the location of the master.

    In three dimensions
    $$
    C =\begin{bmatrix} 
          1 & \widetilde{r} \\
          0 & 1 \\
       \end{bmatrix}.
    $$
    Here $ \widetilde{r}$ is a skew matrix corresponding to the 
    vector $r$, and $0$ and $1$ stand for $3\times3$ zero 
    and identity matrices respectively.

    Further, $\Gamma$ is a diagonal matrix, such that 
    the diagonal entries provide penalty on the difference 
    between the individual degrees of freedom.

    Reference: APPLICATION OF RIGID LINKS  IN STRUCTURAL DESIGN MODELS,
    Sergey Yu. Fialko, International Journal for Computational Civil 
    and Structural Engineering, 13(3) 119-137 (2017). 
    """
    if len(e_x) == 2:
        I = eye(2)
        C = zeros((3, 3))
        C[0:2, 0:2] = I
        C[2, 2] = 1.0
        rx, ry = e_x[0] * h, e_x[1] * h
        C[0, 2] = -ry
        C[1, 2] = +rx
    else:
        I = eye(3)
        C = zeros((6, 6))
        C[0:3, 0:3] = I
        C[3:6, 3:6] = I
        rx, ry, rz = e_x[0] * h, e_x[1] * h, e_x[2] * h
        C[0, 4] = +rz
        C[0, 5] = -ry
        C[1, 3] = -rz
        C[1, 5] = +rx
        C[2, 3] = +ry
        C[2, 4] = -rx
    k = concatenate(
        [
            concatenate([dot(C.T, dot(Gamma, C)), -dot(C.T, Gamma)], axis=1),
            concatenate([-dot(Gamma, C), Gamma], axis=1),
        ],
        axis=0,
    )
    return k


def assemble_stiffness(Kg, member, i, j):
    """
    Assemble rigid link stiffness matrix.

    - `Kg` is the global stiffness matrix,
    - `member` is the rigid link member,
    - `i`, `j` are the joints; the first is the master, the second is the subordinate.
    """
    sect = member["section"]
    Gamma = sect["Gamma"]
    dim = len(i["coordinates"])
    if dim == 2:
        e_x, _, h = geometry.member_2d_geometry(i, j)
    else:
        e_x, _, _, h = geometry.member_3d_geometry(i, j, array([]))
    k = rigid_link_stiffness(e_x, h, Gamma)
    dof = concatenate([i["dof"], j["dof"]])
    return assemble.assemble(Kg, dof, k)
