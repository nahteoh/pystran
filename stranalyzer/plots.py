import matplotlib.pyplot as plt
from stranalyzer.beam import beam_2d_member_geometry, beam_shape_functions, beam_2d_moment
from numpy import linspace, dot, zeros

def plot_setup(m):
    fig = plt.figure()
    if m['dim'] == 3:
        ax = fig.add_subplot(projection='3d')
    else:
        ax = fig.gca()
    return ax

def _plot_members_2d(m):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ci, cj = i['coordinates'], j['coordinates']
        line = plt.plot([ci[0], cj[0]], [ci[1], cj[1]], 'k-')  
    for member in m['beam_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ci, cj = i['coordinates'], j['coordinates']
        line = plt.plot([ci[0], cj[0]], [ci[1], cj[1]], 'k-')  
    return ax
    
def _plot_members_3d(m):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ci, cj = i['coordinates'], j['coordinates']
        line = plt.plot([ci[0], cj[0]], [ci[1], cj[1]], [ci[2], cj[2]], 'k-')
    return ax
    
def plot_members(m):
    if m['dim'] == 3:
        ax = _plot_members_3d(m)
    else:
        ax = _plot_members_2d(m) 
    return ax
    
def _plot_2d_beam(ax, i, j, scale):
    ui, uj = i['displacements'], j['displacements']
    ci, cj = i['coordinates'], j['coordinates']
    e_x, e_y, h = beam_2d_member_geometry(i, j)
    uli = dot(ui[0:2], e_y)
    thi = ui[2]
    ulj = dot(uj[0:2], e_y)
    thj = uj[2]
    n = 20
    xs = zeros(n)
    ys = zeros(n)
    for (s, xi) in enumerate(linspace(-1, +1, n)):
        N = beam_shape_functions(xi, h)
        w = N[0] * uli + (h/2) * N[1] * thi + N[2] * ulj + (h/2) * N[3] * thj
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        xs[s] = x[0] + scale * w * e_y[0]
        ys[s] = x[1] + scale * w * e_y[1]
    ax.plot(xs, ys, 'm-')
        
def plot_deformations(m, scale=1.0):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ui, uj = i['displacements'], j['displacements']
        ci, cj = i['coordinates'], j['coordinates']
        if m['dim'] == 3:
            line = ax.plot([ci[0]+scale*ui[0], cj[0]+scale*uj[0]], [ci[1]+scale*ui[1], cj[1]+scale*uj[1]], [ci[2]+scale*ui[2], cj[2]+scale*uj[2]], 'm-')   
        else:
            line = ax.plot([ci[0]+scale*ui[0], cj[0]+scale*uj[0]], [ci[1]+scale*ui[1], cj[1]+scale*uj[1]], 'm-')   
    for member in m['beam_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        if m['dim'] == 3:
            raise Exception("3D not implemented")  
        else:
            line = _plot_2d_beam(ax, i, j, scale)
    return ax
     
def _plot_member_numbers_2d(m):
    ax = plt.gca()
    for id in m['truss_members'].keys():
        member = m['truss_members'][id]
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ci, cj = i['coordinates'], j['coordinates']
        xm = (ci + cj) / 2.0
        line = ax.text(xm[0], xm[1], str(id))
    return ax
    
def _plot_member_numbers_3d(m):
    ax = plt.gca()
    for id in m['truss_members'].keys():
        member = m['truss_members'][id]
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ci, cj = i['coordinates'], j['coordinates']
        xm = (ci + cj) / 2.0
        line = ax.text(xm[0], xm[1], xm[2], str(id), 'z')
    return ax
    
def plot_member_numbers(m):
    if m['dim'] == 3:
        ax = _plot_member_numbers_3d(m)
    else:
        ax = _plot_member_numbers_2d(m) 
    return ax
     
def _plot_2d_beam_moments(ax, member, i, j, scale):
    e_x, e_y, h = beam_2d_member_geometry(i, j)
    ci, cj = i['coordinates'], j['coordinates']
    n = 7
    xs = zeros(2)
    ys = zeros(2)
    for (s, xi) in enumerate(linspace(-1, +1, n)):
        M = beam_2d_moment(member, i, j, xi)
        x = (1 - xi) / 2 * ci + (1 + xi) / 2 * cj
        #The US convention: moment next to fibers in compression
        xs[0] = x[0]
        xs[1] = x[0] + scale * M * e_y[0]
        ys[0] = x[1]
        ys[1] = x[1] + scale * M * e_y[1]
        ax.plot(xs, ys, 'b-')
    return ax
        
def plot_moments(m, scale=1.0):
    ax = plt.gca()
    for member in m['beam_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        if m['dim'] == 3:
            raise Exception("3D not implemented")  
        else:
            line = _plot_2d_beam_moments(ax, member, i, j, scale)
    return ax
     
def show(m):
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    if m['dim'] == 3:
        ax.set_zlabel('Z')
    plt.show()