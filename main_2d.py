import os
import numpy as np
import fem_2d as fem
import functions_2d as fc
from mesh_process_2d import import_mesh
       
if __name__ == "__main__":
    
    mesh_file = '/home/ana/Downloads/TOPDyn-master/retangulo1x05.IGES' # None
    E = 210e9
    v = 0.3
    rho = 7860
    alpha, beta, eta = 0, 0, 0
    factor = 1000000
    freq = 200
    freq_rsp = [0, 400, 5]

    if mesh_file is not None:
        path = os.path.dirname(os.path.realpath(__file__)) 
        #mesh_file = '/home/ana/Downloads/TOPDyn-master/retangulo1x05.IGES'
        m_file = os.path.join(path, mesh_file)
        coord, connect = import_mesh(m_file)
        ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)
        nelx = len(coord[coord[:, 2] == coord[0, 2]]) - 1
        nely = len(coord[coord[:, 1] == coord[0, 1]]) - 1
        lx = max(coord[:, 1])
        ly = max(coord[:, 2])
    else:
        nelx, nely = 40, 20
        lx, ly = 1, 0.5
        coord, connect, ind_rows, ind_cols = fc.regularmeshQ4(lx, ly, nelx, nely)

    matrix_F = np.array([[1, 0.25, 0, -1, 1]])
    # Create matrix of loads 
    nodes_apply_F = fc.get_nodes_by_coord(coord, matrix_F[:, [0,1]])
    force_matrix = np.empty((nodes_apply_F.shape[0], 4))
    force_matrix[:, 0] = nodes_apply_F
    force_matrix[:, [1,2,3]] = matrix_F[:, [2,3,4]]
    # Create constrain nodes matrix
    restricted_nodes = fc.get_nodes1d(coord, 0, 0.001, 1)
    restri_matrix = np.empty((len(restricted_nodes), 3), dtype='int')
    restri_matrix[:, 0] = restricted_nodes
    restri_matrix[:, [1, 2]] = np.ones((restricted_nodes.shape[0], 2))  
    #
    node_plot = [nodes_apply_F[0], 0, 1]
    save = False

    fem.main(mesh_file, nelx, nely, lx, ly, force_matrix, restri_matrix, E, v, rho, alpha, beta, eta, factor, freq, node_plot=node_plot, freq_rsp=freq_rsp, save=save)





    




    
