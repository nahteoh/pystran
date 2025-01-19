"""

"""

from numpy import array, dot
from numpy.linalg import cross
from context import pystran
from pystran import model
from pystran import section

H = 0.3
B = 0.15
tf = 0.01
tw = 0.005
# p = section.i_beam_points(H, B, tf, tw)
# print(section.area(p))
# print(section.centroid(p))
# p = p - section.centroid(p)
# print(section.centroid(p))
# print(section.inertias(p))

Iy = (B / 12) * H**3 - ((B - tw) / 12) * (H - 2 * tf) ** 3
print(Iy)
Iz = (
    H * B**3 / 12
    - 2 * ((B - tw) / 2) ** 3 * (H - 2 * tf) / 12
    - 2 * ((B - tw) / 2) * (H - 2 * tf) * ((B - tw) / 4 + tw / 2) ** 2
)
print(Iz)
# # Iy = (a3 h / 12) + (b3 / 12) (H - h)

A, Ix, Iy, Iz, J = section.i_beam(H, B, tf, tw)
print(A, Ix, Iy, Iz, J)

print((2 * B * tf**3 + (H - 2 * tf) * tw**3) / 3)
