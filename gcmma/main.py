import numpy as np
import sys
import os
import functions as fc_opt
import plots as plt_opt
sys.path.append(os.getcwd())
import solver_fem_2d.functions_2d as fc
import beam2 as beam 


if __name__ == "__main__":

    nelx, nely = 20, 10
    lx, ly = 0.8, 0.5
    coord, connect, ind_rows, ind_cols = fc.regularmeshQ4(lx, ly, nelx, nely)
    rho = 7860
    E = 210e9
    v = 0.3
    x_min = 0.001
    alpha_par, beta_par, eta_par = 0, 1e-6, 0
    alpha_plot, beta_plot, eta_plot = 0, 1e-6, 0
    p_par = 3
    q_par = 1
    const_func = 100
    matrix_F = np.array([[0.8, 0.25, 0, -1, 10000]])
    # Create matrix of loads 
    nodes_apply_F = fc.get_nodes_by_coord(coord, matrix_F[:, [0,1]])
    force_matrix = np.empty((nodes_apply_F.shape[0], 4))
    force_matrix[:, 0] = nodes_apply_F
    force_matrix[:, [1,2,3]] = matrix_F[:, [2,3,4]]
    # Create constrain nodes matrix
    restri_nodes = fc.get_nodes1d(coord, 0, 0.001, 1)
    restri_matrix = np.empty((len(restri_nodes), 3), dtype='int')
    restri_matrix[:, 0] = restri_nodes
    restri_matrix[:, [1, 2]] = np.ones((restri_nodes.shape[0], 2))
    # Weight at objective function
    n1 = 1
    # Method iterations
    max_iter = 100
    # Factor applied in the radius
    fac_ratio = 2.1
    # If not None is used mode superposition method
    modes = 20
    # Tuple with func_name2 and frequency optimized for func_name2. Associated with weight (1 - n1)
    multiobjective = (None, 0)
    # It can be 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'
    func_name = 'Input Power'
    # Frequency optimized for func_name
    freq1 = 245 + 20
    # Frequency response plot
    freq_rsp = [0, 500, 5]
    # If False use sensitivity filter
    dens_filter = True
    # If True plots the convergence graph for each iteration of the optimization
    each_iter = True
    # Constrain - The first value in the list is used to define the initial value of xval
    constr_func = ['Area']
    constr_values = [50]
    # Plot deformed mesh  
    mesh_deform = True 
    factor = 10000
    # Save plots
    save = False
    # Get code execution time
    timing = False

    beam.main(nelx, nely, lx, ly, func_name, force_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, freq_rsp, dens_filter, each_iter, mesh_deform=mesh_deform, factor=factor, max_iter=max_iter, save=save, timing=timing)
