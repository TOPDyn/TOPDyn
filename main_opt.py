import os
import numpy as np
import functions_opt as fc_opt

if __name__ == "__main__":
    # COM OS EL PASSIVOS
    # If True is used the mma method. If False use the gcmma method.
    mma = True
    
    mesh_file = None
    
    nelx, nely = 200, 100 #0.005 length
    lx, ly = 1, 0.5

    rho = 7860
    E = 210e9
    v = 0.3

    x_min_m = 1e-12
    x_min_k = 1e-9
    alpha_par, beta_par, eta_par = 0, 1e-6, 0 #Não fala nada do eta
    alpha_plot, beta_plot, eta_plot = 0, 1e-6, 0

    p_par = 3
    q_par = 1

    const_func = 100
    
    # Factor applied in the radius
    fac_ratio = 0.021/0.005 #raio/n.el #2.2 #2.1

    # If not None is used mode superposition method
    modes = None

    # Create matrix of loads 
    force_matrix = []
    y_val = np.arange(0.23, 0.27, 0.005)
    for y in y_val:
        aux = [1, y, 0, -1, 1000]
        force_matrix.append(aux)
    
    # Create constrain nodes matrix
    restri_matrix = [[0, 1, 1, 1, 0.001]]
    
    # Weight at objective function
    n1 = 0.8

    # It can be 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio', 'Local Ep', 'Local Ki', 'Local R'
    func_name = 'Input Power' 

    # Frequency optimized for func_name
    freq1 = 50 
    
    # Tuple with func_name2 and frequency optimized for func_name2. Associated with weight (1 - n1)
    multiobjective = ('Compliance', 0)
    
    # Constraint - The first function in the list is used to define the initial value of xval. 'Compliance', 'Local Ep', 'Local Ki', 'Local R' -> (constraint value, frequency)
    constr_func = ['Area', 'Local Ep']
    constr_values = [50, (70, 50)]
    
    passive_coord = ((0.4, 0.6), (0.15, 0.35)) # ((x_initial, x_final), (y_initial, y_final)) or None
   
    # Frequency response plot
    freq_rsp = [1, 1000, 1]
    # If False use sensitivity filter
    dens_filter = True
    # If True plots the convergence graph for each iteration of the optimization
    each_iter = True

    chtol = 1e-4 
    # Plot mesh  
    mesh_deform = False 
    factor = 800
    save = False
    timing = False
    # Method iterations
    max_iter = 171

    nodes = [2]
    number_deltas = 50
    delta_interval = (1e-12, 1e-2) #(low, high)

    #fc_opt.freq_test(mesh_file, nelx, nely, lx, ly, E, v, rho, alpha_par, beta_par, eta_par, 0.5, x_min_m, x_min_k, p_par, q_par, freq_rsp, func_name, const_func, modes, fac_ratio, passive_coord, force_matrix, restri_matrix=restri_matrix)

    #fc_opt.finite_difference(mesh_file, nelx, nely, lx, ly, func_name, force_matrix, restri_matrix, freq1, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, p_par, q_par, passive_coord, nodes, number_deltas, delta_interval)

    fc_opt.exe_opt(mma, mesh_file, nelx, nely, lx, ly, func_name, force_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, passive_coord, freq_rsp, chtol, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)