"""
Define rotation utilities.

"""

from math import pi, cos, sin
import numpy


def rotmat3(rotvec):
    """
    Compute a 3D rotation matrix from a rotation vector.
    """
    R = numpy.zeros((3, 3))
    vn = numpy.linalg.norm(rotvec)
    if vn == 0.0:
        R = numpy.identity(3)
    else:
        rotvecn = (rotvec[0] / vn, rotvec[1] / vn, rotvec[2] / vn)
        ca = cos(vn)
        sa = sin(vn)
        oca = 1.0 - ca
        for j in range(3):
            for i in range(3):
                R[i, j] = oca * rotvecn[i] * rotvecn[j]

        R[0, 1] += sa * -rotvecn[2]
        R[0, 2] += sa * rotvecn[1]
        R[1, 0] += sa * rotvecn[2]
        R[1, 2] += sa * -rotvecn[0]
        R[2, 0] += sa * -rotvecn[1]
        R[2, 1] += sa * rotvecn[0]
        R[0, 0] += ca
        R[1, 1] += ca
        R[2, 2] += ca
    return R
