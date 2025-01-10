import matplotlib.pyplot as plt

def plot_setup(m):
    ax = plt.figure().add_subplot(projection='3d')
    return ax

def plot_members(m):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        xi, xj = i['coordinates'], j['coordinates']
        if m['dim'] == 3:
            line = plt.plot([xi[0], xj[0]], [xi[1], xj[1]], [xi[2], xj[2]], 'k-')
        else:
            line = plt.plot([xi[0], xj[0]], [xi[1], xj[1]], 'k-')   
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    

def plot_deformations(m, scale=1.0):
    ax = plt.gca()
    for member in m['truss_members'].values():
        connectivity = member['connectivity']
        i, j = m['joints'][connectivity[0]], m['joints'][connectivity[1]]
        ui, uj = i['displacements'], j['displacements']
        xi, xj = i['coordinates'], j['coordinates']
        if m['dim'] == 3:
            line = plt.plot([xi[0]+scale*ui[0], xj[0]+scale*uj[0]], [xi[1]+scale*ui[1], xj[1]+scale*uj[1]], [xi[2]+scale*ui[2], xj[2]+scale*uj[2]], 'm-')   
        else:
            line = plt.plot([xi[0]+scale*ui[0], xj[0]+scale*uj[0]], [xi[1]+scale*ui[1], xj[1]+scale*uj[1]], 'm-')   
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
def show():
    plt.show()