"""
Created on Mon Nov 18 08:09:22 2024

@author: pkonl
"""

# import pystran
# from pystran import model
# from pystran import property
# m = model.create()
# model.add_joint(m, 1, [0.1, -0.2])
# print(m['joints'])
# E = 200e9
# A = 0.001
# s1 = property.truss_section('steel', E, A)
# model.add_truss_member(m, 1, [1, 2], s1)
# print(m['truss_members'])

# import pystran
# from pystran import model
# from pystran import property
# m = model.create()
# model.add_joint(m, 1, [0.0, 0.0])
# model.add_support(m['joints'][1], 0, 0.1)
# print(m['joints'][1]['supports'])

# from pystran import model
# from pystran import property
# from pystran import geometry
# m = model.create()
# model.add_joint(m, 1, [0.0, 0.0])
# model.add_joint(m, 2, [1.0, 1.0])
# d = geometry.delt(m['joints'][1]['coordinates'], m['joints'][2]['coordinates'])
# geometry.delt(m['joints'][1]['coordinates'], m['joints'][2]['coordinates'])

# import pystran
# from pystran import model
# from pystran import property
# from pystran import geometry
# m = model.create()
# model.add_joint(m, 1, [0.0, 0.0])
# model.add_joint(m, 2, [1.0, 2.0])
# E = 200e9
# A = 0.001
# s1 = property.truss_section('steel', E, A)
# model.add_truss_member(m, 1, [1, 2], s1)
# i = m['truss_members'][1]['connectivity'][0]
# j = m['truss_members'][1]['connectivity'][1]
# d = geometry.delt(m['joints'][i]['coordinates'], m['joints'][j]['coordinates'])
# print(d)
# L = geometry.len(m['joints'][i]['coordinates'], m['joints'][j]['coordinates'])
# print(L)

# import pystran
# from pystran import model
# from pystran import property
# from pystran import geometry
# from numpy import array, dot, outer
# m = model.create()
# model.add_joint(m, 1, [0.0, 0.0])
# model.add_joint(m, 2, [0.0, -2.0])
# model.add_joint(m, 3, [-1.0, 0.0])
# model.add_joint(m, 4, [+1.0, 0.0])
# model.add_support(m['joints'][1], 0)
# model.add_support(m['joints'][1], 1)
# model.add_support(m['joints'][3], 0)
# model.add_support(m['joints'][3], 1)
# model.add_support(m['joints'][4], 0)
# model.add_support(m['joints'][4], 1)
# E = 200e9
# A = 0.001
# s1 = property.truss_section('steel', E, A)
# model.add_truss_member(m, 1, [1, 2], s1)
# model.add_truss_member(m, 2, [3, 2], s1)
# model.add_truss_member(m, 3, [4, 2], s1)
# model.number_dofs(m)
# print(m['nfreedof'])
# print(m['ntotaldof'])
# print([j['dof'] for j in m['joints'].values()])

# model.solve(m)

import pystran
from pystran import model
from pystran import section
from pystran import geometry
from numpy import array, dot, outer

m = model.create(2)
model.add_joint(m, 1, [0.0, 0.0])
model.add_joint(m, 2, [10.0, 0.0])
model.add_support(m["joints"][1], 0)
model.add_support(m["joints"][1], 1)
model.add_support(m["joints"][1], 2)
E = 200e9
A = 0.001
I = 0.001
s1 = section.beam_2d_section("steel", E, A, I)
model.add_beam_member(m, 1, [1, 2], s1)
model.add_load(m["joints"][2], 2, -1000)
model.number_dofs(m)
print(m["nfreedof"])
print(m["ntotaldof"])
print([j["dof"] for j in m["joints"].values()])

model.solve(m)
