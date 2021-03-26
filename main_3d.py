import numpy as np
import fem_3d as fem
import functions_3d as fc

if __name__ == "__main__":

    nelx, nely, nelz = 10, 10, 10
    lx, ly, lz = 1, 0.5, 0.3

    E = 210e9
    v = 0.3 
    rho = 7860
    alpha, beta, eta = 1e-1, 1e-7, 1e-8
    freq = 380

    coord, connect, ind_rows, ind_cols = fc.regularmeshH8(nelx, nely, nelz, lx, ly, lz)

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

    fem.main(nelx, nely, nelz, lx, ly, lz, load_matrix, restri_matrix, E, v, rho, alpha, beta, eta, factor, freq, freq_rsp, plot_type, complete, amp)
    