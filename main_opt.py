import numpy as np
import functions_2d as fc
import functions_opt as fc_opt
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # If True is used the mma method. If False use the gcmma method.
    mma = True
    nelx, nely = 40, 20
    lx, ly = 1, 0.5
    coord, connect, ind_rows, ind_cols = fc.regularmeshQ4(lx, ly, nelx, nely)
    rho = 7860
    E = 210e9
    v = 0.3
    x_min_m = 1e-6
    x_min_k = 1e-6
    alpha_par, beta_par, eta_par = 0, 1e-5, 0
    alpha_plot, beta_plot, eta_plot = 0, 1e-6, 0
    p_par = 3
    q_par = 1
    const_func = 100
    # Create matrix of loads 
    force_nodes = fc.get_nodes1d(coord, lx, 0.001, 1)
    load_matrix = np.empty((len(force_nodes), 4), dtype='int')
    load_matrix[:, 0] = force_nodes
    load_matrix[:, [1, 2,3]] = [0, -1, 10000]
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
    fac_ratio = 1.1 #2.1
    # If not None is used mode superposition method
    modes = 100
    # Tuple with func_name2 and frequency optimized for func_name2. Associated with weight (1 - n1)
    multiobjective = ('Compliance', 0)
    # It can be 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'
    func_name = 'Input Power'
    # Frequency optimized for func_name
    freq1 = 500
    # Frequency response plot
    freq_rsp = [0, 4000, 50]
    # If False use sensitivity filter
    dens_filter = True
    # If True plots the convergence graph for each iteration of the optimization
    each_iter = True
    # Constraint - The first function in the list is used to define the initial value of xval. 'Compliance' -> (constraint value, frequency)
    constr_func = [ 'Area', 'Compliance']
    constr_values = [50, (20, 0)]
    # Está comentado no while ainda!
    chtol = 1e-4
    # Plot mesh  
    mesh_deform = True 
    factor = 500
    # Save plots
    save = False
    #
    timing = False

    fc_opt.exe_opt(mma, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, freq_rsp, chtol, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)