import os
import numpy as np
import fem_3d as fem

if __name__ == "__main__":
    E = 210e9
    v = 0.3 
    rho = 7860
    alpha, beta, eta = 1e-1, 1e-7, 1e-8
    freq = 380

    path = os.path.dirname(os.path.realpath(__file__)) 
    mesh_file = None
    num_el = None 

    nelx, nely, nelz = 10, 10, 10
    lx, ly, lz = 1, 1, 1
    
    # Constraint nodes in plane x = 0
    restri_matrix =  [{"coord":0, "axis":1, "eps":0.001, "constrain_disp_x":1, "constrain_disp_y":1, "constrain_disp_z":1}]

    # Load 
    load_matrix = [{"x_coord":1, "y_coord":0.3, "z_coord":0.15, "x_direc":0, "y_direc":0, "z_direc":-1, "force":10000}]
    
    # Plot 
    freq_range = [5, 5000, 5]
    factor = 100000
    plot_type = 'deformed'
    complete = True
    amp = 400  

    timing=True

    fem.main(mesh_file, nelx, nely, nelz, lx, ly, lz, load_matrix, restri_matrix, num_el, E, v, rho, alpha, beta, eta, factor, freq, freq_range, plot_type, complete, amp, timing=timing)