import os
import fem_2d as fem
       
if __name__ == "__main__":

    mesh_file = None 

    nelx, nely = 40, 20 
    lx, ly = 1, 0.25

    E = 210e9
    v = 0.3
    rho = 7860

    alpha, beta, eta = 0, 0, 0
    factor = 1000000
    freq = 200
    freq_rsp = [0, 400, 5]

    force_matrix = [[1, 0.25, 0, -1, 1]]

    restri_matrix = [[0, 1, 1, 1, 0.001]]
  
    node_plot = [1, 0.25, 0, 1]
    save = False

    fem.main(mesh_file, nelx, nely, lx, ly, force_matrix, restri_matrix, E, v, rho, alpha, beta, eta, factor, freq, node_plot=node_plot, freq_rsp=freq_rsp, save=save) 
    
    # TESTESS
    #force_matrix = [[1, 0.25, 0, -1, 1], [1, 1, 0, -1, 1], [0.3, 1, 1, 1, 10, 0.001]]

    #restri_matrix = [[0, 2, 0, 1, 0.001], [0, 1, 1, 0, 0.001]] #'valor', coluna (x=1, y=2)
    #restri_matrix = [[0, 2, 0, 1, 0.001], [0, 1, 1, 0, 0.001], [1, 0.5, 1, 1]]
    
    #matrix_F = [[1, 0.25, 0, -1, 1], [1, 1, 0, -1, 1]]

    #matrix_F = [[0, 1, 0.001, 1, 1, 100]]

    #matrix_F = [[0, 1, 0.001, 1, 1, 100], [0, 2, 0.001, 0, 0, 200]]