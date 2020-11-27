from time import time
import functions3D as fc
import numpy as np
import numpy_indexed as npi
import vtkplotter as vt
import pyqtgraph as pg

def all_faces(coord, connect):
    ''' Get vertices of all faces of the mesh '''
    
    nodes_per_face = np.array([connect[:, [1,2,3,4]], connect[:, [5,6,7,8]], \
                                connect[:, [6,7,3,2]], connect[:, [7,8,4,3]], \
                                connect[:, [6,5,1,2]], connect[:, [5,8,4,1]]]).reshape(-1,4)

    ind_faces = npi.indices(coord[:,0], nodes_per_face.flatten()).reshape(-1, 4)
  
    return ind_faces

def free_faces(coord, connect):
    ''' Get vertices of external faces of the mesh '''

    nodes_per_face = np.array([connect[:, [1,2,3,4]], connect[:, [5,6,7,8]], \
                                connect[:, [6,7,3,2]], connect[:, [7,8,4,3]], \
                                connect[:, [6,5,1,2]], connect[:, [5,8,4,1]]]).reshape(-1,4)

    unique, counts = npi.count(nodes_per_face)
    unique = unique[counts<2]  
    ind_faces = npi.indices(coord[:,0], unique.flatten()).reshape(-1, 4) 

    return ind_faces

def build_mesh(verts, ind_faces, scalars=None, timing=False):
    ''' Build the polygonal Mesh '''
    t0 = time()
    mesh = vt.Mesh([verts, ind_faces])
    if scalars is not None:
        mesh.pointColors(scalars, cmap="cool").addScalarBar().normalize()
    else:
        mesh.lineColor('black').lineWidth(1).color('grey').normalize()
    tf = time()
    if timing:
        print("Time to build mesh: " + str(round((tf - t0),6)) + '[s]')

    return mesh

def plot_mesh(mesh, arrows, cones, complete=True, mesh2=None):
    
    if mesh2 is not None:
        if complete:
            vt.show(mesh, arrows, cones, mesh2, __doc__, viewup='z', axes=4)
        else:
            vt.show(mesh, mesh2, __doc__, viewup='z', axes=4)
    else:  
        if complete:
            vt.show(mesh, arrows, cones, __doc__, viewup='z', axes=4)
        else:
            vt.show(mesh, __doc__, viewup='z', axes=4)

def animation(coord, connect, amp, disp_vector):
    '''
    amp: amplitude
    '''
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

def build_arrows(load_matrix, normalized_coord, timing=False):
    t0 = time()
    force_coord = get_normalized_coord(load_matrix[:, 0], normalized_coord)
    delta = 0.2 * normalized_coord.max()
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

def build_cones(restricted, normalized_coord, timing=False):
    t0 = time()
    restri_coord = get_normalized_coord(restricted[:, 0], normalized_coord)
    max_value = normalized_coord.max()
    ratio = 0.02 * max_value
    height = 0.05 * max_value
    delta = height/2
    cones = []

    for i in range(restricted.shape[0]):
        if restricted[i, 1] == 1:
            cone = vt.Cone([restri_coord[i, 0] - delta, restri_coord[i, 1], restri_coord[i, 2]], r=ratio, height=height, axis=[1, 0, 0], c='red', alpha=1, res=48)
            cones.append(cone)

        if restricted[i, 2] == 1:
            cone = vt.Cone([restri_coord[i, 0], restri_coord[i, 1] - delta, restri_coord[i, 2]], r=ratio, height=height, axis=[0, 1, 0], c='green', alpha=1, res=48)
            cones.append(cone)

        if restricted[i, 3] == 1:
            cone = vt.Cone([restri_coord[i, 0], restri_coord[i, 1], restri_coord[i, 2] - delta], r=ratio, height=height, axis=[0, 0, 1], c='blue', alpha=1, res=48)
            cones.append(cone)
    tf = time()
    if timing:
        print("Time to build cones: " + str(round((tf - t0),6)) + '[s]')

    return cones

def get_normalized_coord(node, normalized_coord):

    idx_node = (node - 1).astype('int')
    coord = normalized_coord[idx_node, :]
    return coord

def plot_freqresponse(freq_range, delta, disp_vector):

    win = pg.GraphicsLayoutWidget(show=True, title="MMA")
    win.resize(1000,600)
    #p = win.addPlot(title="Convergence")
    p = win.addPlot()
    p.addLegend(labelTextColor=(0,0,0), offset=(800,10))
    x = np.arange(freq_range[0], freq_range[1] + 1, delta)
    y = 10 * np.log10(abs(disp_vector))
    p.plot(x, y, width=80, pen=(0,0,0))
    colors = [(43,174,179), (64,66,114), (255,110,60), (255,215,75)]
    p.setLabel('left', 'displacement, U [N]')
    p.setLabel('bottom', "frequency, f [Hz]") 

    return win

def select_nodes(coord, connect, elements):
    ''' Select nodes to build the mesh '''
    elements = np.array(elements) - 1
    nodes = connect[elements, 1:].flatten()
    index = npi.indices(coord[:, 0], np.unique(nodes))
    
    return coord[index, :], connect[elements,:]