import matplotlib.pyplot as plt

def plot_setup(m):
    if m['dim'] == 3:
        ax = plt.figure().add_subplot(projection='3d')
    else:
        ax = plt.gca()
    return ax

def _plot_members_2d(m):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        xi, xj = i['coordinates'], j['coordinates']
        line = plt.plot([xi[0], xj[0]], [xi[1], xj[1]], 'k-')   
    return ax
    
def _plot_members_3d(m):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        xi, xj = i['coordinates'], j['coordinates']
        line = plt.plot([xi[0], xj[0]], [xi[1], xj[1]], [xi[2], xj[2]], 'k-')
    return ax
    
def plot_members(m):
    if m['dim'] == 3:
        ax = _plot_members_3d(m)
    else:
        ax = _plot_members_2d(m) 
    return ax
    

def plot_deformations(m, scale=1.0):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ui, uj = i['displacements'], j['displacements']
        xi, xj = i['coordinates'], j['coordinates']
        if m['dim'] == 3:
            line = ax.plot([xi[0]+scale*ui[0], xj[0]+scale*uj[0]], [xi[1]+scale*ui[1], xj[1]+scale*uj[1]], [xi[2]+scale*ui[2], xj[2]+scale*uj[2]], 'm-')   
        else:
            line = ax.plot([xi[0]+scale*ui[0], xj[0]+scale*uj[0]], [xi[1]+scale*ui[1], xj[1]+scale*uj[1]], 'm-')   
    return ax
     
def _plot_member_numbers_2d(m):
    ax = plt.gca()
    for id in m['truss_members'].keys():
        member = m['truss_members'][id]
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        xi, xj = i['coordinates'], j['coordinates']
        xm = (xi + xj) / 2.0
        line = ax.text(xm[0], xm[1], str(id))
    return ax
    
def _plot_member_numbers_3d(m):
    ax = plt.gca()
    for id in m['truss_members'].keys():
        member = m['truss_members'][id]
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        xi, xj = i['coordinates'], j['coordinates']
        xm = (xi + xj) / 2.0
        line = ax.text(xm[0], xm[1], xm[2], str(id), 'z')
    return ax
    
def plot_member_numbers(m):
    if m['dim'] == 3:
        ax = _plot_member_numbers_3d(m)
    else:
        ax = _plot_member_numbers_2d(m) 
    return ax
    return ax
    
def show(m):
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    if m['dim'] == 3:
        ax.set_zlabel('Z')
    plt.show()