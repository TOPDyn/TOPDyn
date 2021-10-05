from time import time
import numpy as np
import pyqtgraph as pg
import vedo as vt
import functions_3d as fc
import numpy_indexed as npi

def all_faces(coord, connect):
    """ Gets vertices of all faces of the mesh.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
    
    Returns:
        Corresponding nodes.
    """   
    nodes_per_face = np.array([connect[:, [1,2,3,4]], connect[:, [5,6,7,8]], \
                                connect[:, [6,7,3,2]], connect[:, [7,8,4,3]], \
                                connect[:, [6,5,1,2]], connect[:, [5,8,4,1]]]).reshape(-1,4)

    ind_faces = npi.indices(coord[:,0], nodes_per_face.flatten()).reshape(-1, 4)
    return ind_faces

def free_faces(coord, connect):
    """ Gets vertices of external faces of the mesh.
    
    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
    
    Returns:
        Corresponding nodes.
    """
    nodes_per_face = np.array([connect[:, [1,2,3,4]], connect[:, [5,6,7,8]], \
                                connect[:, [6,7,3,2]], connect[:, [7,8,4,3]], \
                                connect[:, [6,5,1,2]], connect[:, [5,8,4,1]]]).reshape(-1,4)

    unique, counts = npi.count(nodes_per_face)
    unique = unique[counts<2]  
    ind_faces = npi.indices(coord[:,0], unique.flatten()).reshape(-1, 4) 
    return ind_faces

def build_mesh(verts, ind_faces, scalars=None, timing=False):
    """ Builds the polygonal Mesh.
    
    Args:
        verts (:obj:`numpy.array`): Vertices of the elements.
        ind_faces (:obj:`numpy.array`): Connectivity of the faces of the elements.
        scalars (:obj:`int`, optional): Values to add the scalar bar. Defaults to None.
        timing (:obj:`bool`, optional): If True shows the time to build the mesh. Defaults to False.
    
    Returns:
        Build an instance of the Mesh object from the Vedo library.
    """
    t0 = time()
    mesh = vt.Mesh([verts, ind_faces])
    if scalars is not None:
        mesh.pointColors(scalars, cmap="cool").addScalarBar()#.normalize()
    else:
        mesh.lineColor('black').lineWidth(1).color('grey')#.normalize()
    tf = time()
    if timing:
        print("Time to build mesh: " + str(round((tf - t0),6)) + '[s]')
    return mesh

def plot_mesh(mesh, arrows, cones, complete=True, mesh2=None):
    """ Plot mesh.

    Args:
        mesh (vedo.pointcloud.Points): Mesh object.
        arrows (vedo.shapes.Arrows): Load arrows.
        cones (vedo.shapes.Cones): Cones to represent constrain nodes.
        complete (:obj:`bool`, optional): If true plot mesh with loads and constrain nodes. Defaults to True.
        mesh2 (:obj:`vedo.pointcloud.Points`, optional): Mesh without faces. Defaults to None.
    """
    vp = vt.Plotter()

    vp += mesh
    if complete:
        vp += arrows
        vp += cones
    if mesh2 is not None:
        vp += mesh2
    return vp
    
def animation(coord, connect, amp, disp_vector):
    """ Plot deformed mesh animation.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        amp (:obj:`int`): Amplitude. 
        disp_vector (:obj:`numpy.array`): Displacement.
    """
    vt.printc("Press F1 to exit.", c="red", invert=1)
    vp = vt.Plotter(axes=0, interactive=0)
    vp += __doc__
    ind_faces_U = free_faces(coord, connect)
    aux = np.linspace(-1, 1, num=400)
    factor = amp * np.sin(aux * 2 * np.pi)
    while True:
        for f in factor:
            coord3D_U = fc.apply_U(disp_vector, coord, factor=f)
            verts_U = coord3D_U[:, 1:]
            mesh = build_mesh(verts_U, ind_faces_U)
            vp.clear()
            vp.add(mesh)
            vp.show()

def build_arrows(load_matrix, normal_coord, timing=False):
    """ Builds load arrows.
        
    Args:  
        load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value.   
        normal_coord (:obj:`numpy.array`): mesh coordinates.
        timing (:obj:`bool`, optional): If True shows the time to build the load arrows. Defaults to False.
       
    Returns:
        List with the load arrows.
    """
    t0 = time()
    force_coord = get_normalized_coord(load_matrix[:, 0], normal_coord)
    delta = 0.2 * normal_coord.max()
    arrows_start, arrows_end = [], []

    for i in range(load_matrix.shape[0]):
        
        if load_matrix[i, 1] == 1:
            arrows_in = list(force_coord[i, :])
            arrows_out = [force_coord[i, 0] + delta, force_coord[i, 1], force_coord[i, 2]]
            arrows_start.append(arrows_in)
            arrows_end.append(arrows_out)

        if load_matrix[i, 2] == 1:
            arrows_in = list(force_coord[i, :])
            arrows_out = [force_coord[i, 0], force_coord[i, 1] + delta, force_coord[i, 2]]
            arrows_start.append(arrows_in)
            arrows_end.append(arrows_out)

        if load_matrix[i, 3] == 1:
            arrows_in = list(force_coord[i, :])
            arrows_out = [force_coord[i, 0], force_coord[i, 1], force_coord[i, 2] - delta]
            arrows_start.append(arrows_in)
            arrows_end.append(arrows_out)
    
        if load_matrix[i, 1] == -1:
            arrows_in = [force_coord[i, 0] + delta, force_coord[i, 1], force_coord[i, 2]]
            arrows_out = list(force_coord[i, :])
            arrows_start.append(arrows_in)
            arrows_end.append(arrows_out)

        if load_matrix[i, 2] == -1:
            arrows_in = [force_coord[i, 0], force_coord[i, 1] + delta, force_coord[i, 2]]
            arrows_out = list(force_coord[i, :])
            arrows_start.append(arrows_in)
            arrows_end.append(arrows_out)

        if load_matrix[i, 3] == -1:
            arrows_in = [force_coord[i, 0], force_coord[i, 1], force_coord[i, 2] - delta]
            arrows_out = list(force_coord[i, :])
            arrows_start.append(arrows_in)
            arrows_end.append(arrows_out)

    all_arrows = vt.Arrows(arrows_start, arrows_end, c='darkorange')
    tf = time()
    if timing:
        print("Time to build arrows: " + str(round((tf - t0),6)) + '[s]')
    return all_arrows

def build_cones(restri_nodes, normal_coord, timing=False):
    """ Builds cones to indicate constrain nodes.
        
    Args:  
        restri_nodes (:obj:`numpy.array`): Constrain nodes.  
        normal_coord (:obj:`numpy.array`): mesh coordinates.
        timing (:obj:`bool`, optional): If True shows the time to build the load arrows. Defaults to False.
       
    Returns:
        List with the cones.
    """
    t0 = time()
    restri_coord = get_normalized_coord(restri_nodes[:, 0], normal_coord)
    max_value = normal_coord.max()
    ratio = 0.02 * max_value
    height = 0.05 * max_value
    delta = height/2
    cones = []

    for i in range(restri_nodes.shape[0]):
        if restri_nodes[i, 1] == 1:
            cone = vt.Cone([restri_coord[i, 0] - delta, restri_coord[i, 1], restri_coord[i, 2]], r=ratio, height=height, axis=[1, 0, 0], c='red', alpha=1, res=48)
            cones.append(cone)

        if restri_nodes[i, 2] == 1:
            cone = vt.Cone([restri_coord[i, 0], restri_coord[i, 1] - delta, restri_coord[i, 2]], r=ratio, height=height, axis=[0, 1, 0], c='green', alpha=1, res=48)
            cones.append(cone)

        if restri_nodes[i, 3] == 1:
            cone = vt.Cone([restri_coord[i, 0], restri_coord[i, 1], restri_coord[i, 2] - delta], r=ratio, height=height, axis=[0, 0, 1], c='blue', alpha=1, res=48)
            cones.append(cone)
    tf = time()
    if timing:
        print("Time to build cones: " + str(round((tf - t0),6)) + '[s]')
    return cones

def get_normalized_coord(nodes, normal_coord):
    """ Normalizes coordinates.

    Args:
        nodes (:obj:`numpy.array`): External face nodes.
        normal_coord (:obj:`numpy.array`): Normalized coordinate.

    Returns:
       Normalized coordinate of the external face nodes.
    """
    idx_node = (nodes - 1).astype('int')
    coord = normal_coord[idx_node, :]
    return coord

def plot_freqresponse(freq_range, disp_vector):
    """ Plot the frequency response.
            
    Args:    
        freq_range (:obj:`list`): Range of frequencies analyzed.
            
            * First value is the minimum frequency.
            * Second value is the maximum frequency.
            * Third value is the step between each calculation of the function. 
        disp_vector (:obj:`numpy.array`): Displacement.
       
    Returns:
        Window with the graph.
    """
    win = pg.GraphicsLayoutWidget(show=True, title="MMA")
    win.resize(1000,600)
    #p = win.addPlot(title="Convergence")
    p = win.addPlot()
    p.addLegend(labelTextColor=(0,0,0), offset=(800,10))
    x = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
    y = 10 * np.log10(abs(disp_vector))
    p.plot(x, y, width=80, pen=(0,0,0))
    colors = [(43,174,179), (64,66,114), (255,110,60), (255,215,75)]
    p.setLabel('left', 'displacement, U [N]')
    p.setLabel('bottom', "frequency, f [Hz]") 
    return win, p

def select_nodes(coord, connect, nodes):
    """ Selects unique nodes to build the mesh.
    
    Args:
        coord (:obj:`numpy.array`): mesh coordinates.
        connect (:obj:`numpy.array`): Element connectivity.
        nodes (:obj:`numpy.array`): Nodes to select.

    Returns:
        A tuple with coordinates and connectivity for the selected nodes.
    """
    nodes = np.array(nodes) - 1
    nodes = connect[nodes, 1:].flatten()
    index = npi.indices(coord[:, 0], np.unique(nodes))
    return coord[index, :], connect[nodes,:]