import numpy as np
import FEM2D as fem
import functions2D as fc
       
if __name__ == "__main__":
  
    nelx, nely = 20, 10
    lx, ly = 1, 0.5
    E = 210e9
    v = 0.3
    rho = 7860
    alpha, beta, eta = 0, 0, 0
    factor = 1 
    freq = 200
    freq_rsp = [400, 5]
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

    fem.main(nelx, nely, lx, ly, force_matrix, None, E, v, rho, alpha, beta, eta, factor, freq, freq_rsp=freq_rsp)





    




    
