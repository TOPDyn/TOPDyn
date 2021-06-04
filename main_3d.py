import os
import numpy as np
import fem_3d as fem
import functions_3d as fc
from mesh_process_3d import import_mesh

if __name__ == "__main__":

    E = 210e9
    v = 0.3 
    rho = 7860
    alpha, beta, eta = 1e-1, 1e-7, 1e-8
    freq = 380

    path = os.path.dirname(os.path.realpath(__file__)) 
    mesh_file = '/home/ana/Downloads/TOPDyn-master/z_cubo1x1x1.IGES'
    num_el = 10 # É o número de elementos em uma direção, não o total. Default é 30 igual ao do gmsh
    if mesh_file is not None:
        m_file = os.path.join(path, mesh_file)
        coord, connect = import_mesh(m_file, num_el)
        lx = max(coord[:, 1])
        ly = max(coord[:, 2])
        lz = max(coord[:, 3])
        nelx = len(coord[np.logical_and(coord[:, 2] == coord[0, 2], coord[:, 3] == coord[0, 3])]) - 1
        nely = len(coord[np.logical_and(coord[:, 1] == coord[0, 1], coord[:, 3] == coord[0, 3])]) - 1
        nelz = len(coord[np.logical_and(coord[:, 1] == coord[0, 1], coord[:, 2] == coord[0, 2])]) - 1
    else:
        nelx, nely, nelz = 10, 10, 10
        lx, ly, lz = 1, 1, 1
        coord, connect, _, _ = fc.regularmeshH8(nelx, nely, nelz, lx, ly, lz)

    # Restricted nodes in plane x = 0
    restri_matrix = np.empty(((nely + 1) * (nelz + 1), 4))
    mask = coord[:, 1] == 0.
    restri_matrix[:, 0] = (coord[mask, 0]).astype('int')
    restri_matrix[:, [1,2,3]] = np.ones(((nely + 1) * (nelz + 1), 3), dtype='int')
    # Force applied in coordinate (1, 0.3, 0.15)
    load_matrix = np.empty((1,5))
    load_matrix[0, 0] = fc.get_nodes_by_coord(coord, np.array([[1, 0.3, 0.15]]))
    load_matrix[0, [1,2,3,4]] = np.array([[0, 0, -1, 10000]])
    # Plot 
    freq_rsp = [0, 400, 5]
    factor = 100000
    plot_type = 'deformed'
    complete = True
    amp = 400  

    fem.main(mesh_file, num_el, nelx, nely, nelz, lx, ly, lz, load_matrix, restri_matrix, E, v, rho, alpha, beta, eta, factor, freq, freq_rsp, plot_type, complete, amp)
    