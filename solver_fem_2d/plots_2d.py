from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as cl
import matplotlib.cm as cm
import matplotlib.tri as tri

def show_nodes(coord):
    """ Plot nodes of mesh.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
    """
    ax = plt.axes()
    ax.scatter(coord[:,1], coord[:,2])

    for i in range(coord.shape[0]):
        ax.annotate(coord[i,0], (coord[i,1], coord[i,2]))

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

def build_collection(coord, connect, disp_vector=None, timing=False):
    """ Build quad mesh.
    
    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        disp_vector (:obj:`numpy.array`): Displacement.
        timing (:obj:`bool`, optional): If True shows the time to build the mesh. Defaults to False.
     
    Returns:
        Matplotlib collection object.
    """
    t0 = time()
    x, y = coord[:,1], coord[:,2]
    xy = np.c_[x,y]
    squares = np.asarray(connect)
    verts = xy[squares - 1]
    pc = cl.PolyCollection(verts)
    if disp_vector is not None:
        pc.set_array(disp_vector)
        pc.set_cmap("viridis")
    else:
        pc.set_edgecolor("black")
        pc.set_facecolor("None")
    tf = time()
    if timing:
        print("Time to build collection: " + str(round((tf - t0),6)) + '[s]')

    return pc

def plot_collection(lx, ly, coord, pc, load_matrix=None, restri_matrix=None, save=False):
    """ Plot mesh, force arrows and constrain nodes. 

    Args:
        lx (:obj:`int`): X-axis length.
        ly (:obj:`int`): Y-axis length.
        coord (:obj:`numpy.array`): Coordinates of the element.
        pc (matplotlib.collections): Matplotlib collection object.
        load_matrix (:obj:`numpy.array`, optional): The columns are respectively node, x direction, y direction, force value. 
        restri_matrix (:obj:`numpy.array`, optional)= The columns are respectively node, x direction, y direction. Defaults to None. 
        save (:obj:`bool`, optional): True for save the graphic in PNG. Defaults to False.

    Returns:
        Matplotlib axes object.
    """
    fig, ax = plt.subplots()
    ax.add_collection(pc)
    ax.autoscale()
    ax.set_aspect('equal')
    x, y = coord[:, 1], coord[:, 2]
    ax.plot(x, y, ls = "", color = "black")
    max_size, min_size = get_size(lx, ly)
    size = 0.06 * coord[:, max_size].max() - coord[:, max_size].min()
    # Plot force arrows
    if load_matrix is not None:
        ax = plot_force(ax, coord, load_matrix, min_size, size)
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
    """ Get columns with maximum and minimum length.

    Args:
        lx (:obj:`int`): X-axis length.
        ly (:obj:`int`): Y-axis length.

    Returns:
        Column indexes with maximum and minimum length, respectively.
    """
    if lx > ly:
        max_size = 1
        min_size = 2
    else:
        max_size = 2
        min_size = 1

    return max_size, min_size

def plot_force(ax, coord, load_matrix, min_size, size):
    """ Add load vector to the plot of deformed mesh.
        
    Args:    
        ax (matplotlib.axes.Axes): Matplotlib axes object.
        coord (:obj:`numpy.array`): mesh coordinates.
        load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value. 
        min_size (:obj:`float`): Value to define the width of load arrow.
        size (:obj:`float`): Load arrow length.
    
    Returns:
        Matplotlib axes object.
    """
    factor = coord[:, min_size].max() - coord[:, min_size].min()
    ind = (load_matrix[:, 0] - 1).astype('int')
    for i in range(load_matrix.shape[0]):
        # 
        if load_matrix[i, 1] == 1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
        if load_matrix[i, 1] == -1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], -size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
        # 
        if load_matrix[i, 2] == -1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], 0, -size, shape='full', length_includes_head=True, width = factor * 0.01)
        #
        if load_matrix[i, 2] == 1:
            ax.arrow(coord[ind[i], 1], coord[ind[i], 2], 0, size, shape='full', length_includes_head=True, width = factor * 0.01)
    
    return ax

def plot_restr_nodes(ax, coord, restri_matrix):
    """ Add contrain nodes to the plot of deformed mesh.
        
    Args:    
        ax (matplotlib.axes.Axes): Matplotlib axes object.
        coord (:obj:`numpy.array`): mesh coordinates.
        restri_matrix (numpy.array)= The columns are respectively node, x direction, y direction.
       
    Returns:
        Matplotlib axes object.
    """
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

def plot_freqresponse(freq_range, delta, disp_vector, save=False):
    """ Plot the frequency response.
            
    Args:    
        freq_range (:obj:`list`): Range of frequencies analyzed.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (:obj:`int`): Step between each calculation of the function. 
        disp_vector (:obj:`numpy.array`): Displacement.
        save (:obj:`bool`, optional): True for save the graphic in PNG. Defaults to False.
       
    Returns:
        Matplotlib axes object.
    """
    fig, ax = plt.subplots()
    x = np.arange(freq_range[0], freq_range[1] + 1, delta)
    y = 10 * np.log10(abs(disp_vector))
    ax.plot(x, y)
    ax.set_xlabel('frequency [Hz]', fontsize=16)
    ax.set_ylabel('displacement [N]', fontsize=16)
    ax.set_xlim(0, x[-1])
    ax.grid()
    if save:
        plt.savefig("freq_range.png")

    return ax

def many_disp(freq_range, delta, disp_vector1, disp_vector2, disp_vector3):
    """ Plot the frequency response of three displacement vectors.
            
    Args:    
        freq_range (:obj:`list`): Range of frequencies analyzed.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (:obj:`int`): Step between each calculation of the function. 
        disp_vector1 (:obj:`numpy.array`): Displacement 1.
        disp_vector2 (:obj:`numpy.array`): Displacement 2.
        disp_vector3 (:obj:`numpy.array`): Displacement 3.
    """
    x = np.arange(freq_range[0], freq_range[1] + 1, delta)
    y1 = 10 * np.log10(abs(disp_vector1))
    y2 = 10 * np.log10(abs(disp_vector2))
    y3 = 10 * np.log10(abs(disp_vector3))

    plt.plot(x, y1, 'r')
    plt.plot(x, y2, 'g')
    plt.plot(x, y3, 'b')

    plt.xlabel('frequency [Hz]', fontsize=16)
    plt.ylabel('displacement [N]', fontsize=16)
    plt.xlim(0, x[-1])
    plt.grid()
    plt.legend(['0', '0.5', '1'])
    plt.show()
