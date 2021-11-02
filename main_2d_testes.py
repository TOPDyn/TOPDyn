import os
import fem_2d as fem
       
if __name__ == "__main__":

    mesh_file = None 

    nelx, nely = 10, 5 
    lx, ly = 1, 0.5

    E = 210e9
    v = 0.3
    rho = 7860

    alpha, beta, eta = 0, 0, 0
    factor = 1000000
    freq = 200
    freq_rsp = [2, 400, 2]

    force_matrix = [{"x_coord":1, "y_coord":0.25, "x_direc":0, "y_direc":-1, "force":1}]

    restri_matrix = [{"coord":0, "axis":1, "eps":0.001, "constrain_disp_x":1, "constrain_disp_y":1}]
  
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