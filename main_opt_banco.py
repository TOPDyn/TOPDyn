import numpy as np
import functions_2d as fc
import functions_opt as fc_opt

if __name__ == "__main__":
    # If True is used the mma method. If False use the gcmma method.
    mma = True
    nelx, nely = 5, 10
    lx, ly = 0.5, 1
    rho = 7860
    E = 210e9
    v = 0.3
    x_min_m = 1e-12
    x_min_k = 1e-9
    alpha_par, beta_par, eta_par = 0, 1e-5, 0
    alpha_plot, beta_plot, eta_plot = 0, 1e-6, 0
    p_par = 3
    q_par = 1
    const_func = 100

    # Create matrix of loads 
    load_matrix = [{"coord":ly, "axis":2, "x_direc":0, "y_direc":-1, "force":10000, "eps":0.001}]
    
    # Create constrain nodes matrix
    restri_matrix = [{"coord":0, "axis":2, "eps":0.001, "constrain_disp_x":1, "constrain_disp_y":1}]

    # Weight at objective function
    n1 = 1
    # Method iterations
    max_iter = 3
    # Factor applied in the radius
    fac_ratio = 2.2 #2.1
    # If not None is used mode superposition method
    modes = None
    # Tuple with func_name2 and frequency optimized for func_name2. Associated with weight (1 - n1)
    multiobjective = ('compliance', 0)
    # It can be 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'
    func_name = 'input_power'
    # Frequency optimized for func_name
    freq1 = 10
    # Frequency response plot
    freq_rsp = [5, 500, 5]
    # If False use sensitivity filter
    dens_filter = True
    # If True plots the convergence graph for each iteration of the optimization
    each_iter = True
    # Constraint - The first function in the list is used to define the initial value of xval. 'Compliance' -> (constraint value, frequency)
    constr_func = [ 'area']
    constr_values = [30]
    
    passive_coord = ((0, 0.5), (0.95, 1)) # ((x_initial, x_final), (y_initial, y_final)) or None

    # Plot mesh  
    mesh_deform = True 
    factor = 0
    # Save plots
    save = True
    timing = False

    fc_opt.exe_opt(mma, None, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, passive_coord, freq_rsp, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)