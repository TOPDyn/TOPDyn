from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as cl
import matplotlib.cm as cm
import matplotlib.tri as tri

def show_nodes(coord):
    ''' Plot nodes of mesh'''
    ax = plt.axes()
    ax.scatter(coord[:,1], coord[:,2])

    for i in range(coord.shape[0]):
        ax.annotate(coord[i,0], (coord[i,1], coord[i,2]))

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

def build_collection(coord, connect, displacement=None, timing=False):
    ''' Build quad mesh'''
    t0 = time()
    x, y = coord[:,1], coord[:,2]
    xy = np.c_[x,y]
    squares = np.asarray(connect)
    verts = xy[squares - 1]
    pc = cl.PolyCollection(verts)
    if displacement is not None:
        pc.set_array(displacement)
        pc.set_cmap("viridis")
    else:
        pc.set_edgecolor("black")
        pc.set_facecolor("None")
    tf = time()
    if timing:
        print("Time to build mesh: " + str(round((tf - t0),6)) + '[s]')

    return pc

def plot_collection(lx, ly, coord, pc, force_matrix=None, restri_matrix=None, save=False):
    ''' Plot mesh, force arrows and restricted nodes '''
    fig, ax = plt.subplots()
    ax.add_collection(pc)
    ax.autoscale()
    ax.set_aspect('equal')
    x, y = coord[:, 1], coord[:, 2]
    ax.plot(x, y, ls = "", color = "black")
    max_size, min_size = get_size(lx, ly)
    size = 0.06 * coord[:, max_size].max() - coord[:, max_size].min()
    # Plot force arrows
    if force_matrix is not None:
        ax = plot_force(ax, coord, force_matrix, min_size, size)
    # Plot restrict nodes 
    if restri_matrix is not None:
        ax = plot_restr_nodes(ax, coord, restri_matrix)
    #
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_xlim(coord[:, 1].min() - size, coord[:, 1].max() + size)
    ax.set_ylim(coord[:, 2].min() - size, coord[:, 2].max() + size)
    if save:
        plt.savefig("mesh.png")
    return ax

def get_size(lx, ly):
    if lx > ly:
        max_size = 1
        min_size = 2
    else:
        max_size = 2
        min_size = 1

    return max_size, min_size

def plot_force(ax, coord, force_matrix, min_size, size):
    factor = coord[:, min_size].max() - coord[:, min_size].min()
    ind = (force_matrix[:, 0] - 1).astype('int')
    for i in range(force_matrix.shape[0]):
        # 
        if force_matrix[i, 1] == 1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
        if force_matrix[i, 1] == -1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], -size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
        # 
        if force_matrix[i, 2] == -1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], 0, -size, shape='full', length_includes_head=True, width = factor * 0.01)
        #
        if force_matrix[i, 2] == 1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], 0, size, shape='full', length_includes_head=True, width = factor * 0.01)
    return ax

def plot_restr_nodes(ax, coord, restri_matrix):
    both_restri = (restri_matrix[:, 1] == 1) & (restri_matrix[:, 2] == 1)
    ind = restri_matrix[both_restri, 0] - 1
    ax.scatter(coord[ind, 1], coord[ind, 2], marker=(3, 0, 0), s=120, color = 'red', linestyle='None')
    dist = 0.01 * (coord[1,0] - coord[0,0])
    ax.scatter(coord[ind, 1] - dist, coord[ind, 2], marker=(3, 0, 270), s=120, color = 'green', linestyle='None')
    
    x_restri = (restri_matrix[:, 1] == 1) & (restri_matrix[:, 2] == 0)
    ind = restri_matrix[x_restri, 0] - 1
    ax.scatter(coord[ind, 1], coord[ind, 2], marker=(3, 0, 270), s=120, color = 'green', linestyle='None')
    
    y_restri = (restri_matrix[:, 1] == 0) & (restri_matrix[:, 2] == 1)
    ind = restri_matrix[y_restri, 0] - 1
    plt.scatter(coord[ind, 1], coord[ind, 2], marker=(3, 0, 0), s=120, color = 'red', linestyle='None')

    return ax

def plot_freqresponse(freq_range, delta, vector, save=False):
    fig, ax = plt.subplots()
    x = np.arange(freq_range[0], freq_range[1] + 1, delta)
    y = 10 * np.log10(abs(vector))
    ax.plot(x, y)
    ax.set_xlabel('frequency [Hz]', fontsize=16)
    ax.set_ylabel('displacement [N]', fontsize=16)
    ax.set_xlim(0, x[-1])
    ax.grid()
    if save:
        plt.savefig("freq_range.png")
    return ax

def many_frequencies(freq_range, delta, vector1, vector2, vector3):
    x = np.arange(freq_range[0], freq_range[1] + 1, delta)
    y1 = 10 * np.log10(abs(vector1))
    y2 = 10 * np.log10(abs(vector2))
    y3 = 10 * np.log10(abs(vector3))

    plt.plot(x, y1, 'r')
    plt.plot(x, y2, 'g')
    plt.plot(x, y3, 'b')

    plt.xlabel('frequency [Hz]', fontsize=16)
    plt.ylabel('displacement [N]', fontsize=16)
    plt.xlim(0, x[-1])
    plt.grid()
    plt.legend(['0', '0.5', '1'])
    plt.show()
