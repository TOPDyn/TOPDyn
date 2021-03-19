import numpy as np
import sys
import os
sys.path.append(os.getcwd())
import solver_fem_2d.functions_2d as fc
import beam2 as beam 
import functions as fc_opt
import matplotlib.pyplot as plt

if __name__ == "__main__":

    nelx, nely = 20, 10
    lx, ly = 3, 0.5
    coord, connect, ind_rows, ind_cols = fc.regularmeshQ4(lx, ly, nelx, nely)
    rho = 7860
    E = 210e9
    v = 0.3
    x_min = 1e-6
    alpha_par, beta_par, eta_par = 0, 1e-6, 0
    alpha_plot, beta_plot, eta_plot = 0, 1e-6, 0
    p_par = 3
    q_par = 1
    const_func = 100
    # Create matrix of loads 
    force_nodes = fc.get_nodes1d(coord, 3, 0.001, 1)
    force_matrix = np.empty((len(force_nodes), 4), dtype='int')
    force_matrix[:, 0] = force_nodes
    force_matrix[:, [1, 2,3]] = [-1, 0, 10000]
    # Create constrain nodes matrix
    restri_nodes = fc.get_nodes1d(coord, 0, 0.001, 1)
    restri_matrix = np.empty((len(restri_nodes), 3), dtype='int')
    restri_matrix[:, 0] = restri_nodes
    restri_matrix[:, [1, 2]] = np.ones((restri_nodes.shape[0], 2))
    # Weight at objective function
    n1 = 0.9
    # Method iterations
    max_iter = 100
    # Factor applied in the radius
    fac_ratio = 2.1
    # If not None is used mode superposition method
    modes = 400
    # Tuple with func_name2 and frequency optimized for func_name2. Associated with weight (1 - n1)
    multiobjective = ('Compliance', 0)
    # It can be 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'
    func_name = 'Input Power'
    # Frequency optimized for func_name
    freq1 = 150 #215
    # Frequency response plot
    freq_rsp = [0, 4000, 50]
    # If False use sensitivity filter
    dens_filter = True
    # If True plots the convergence graph for each iteration of the optimization
    each_iter = True
    # Constraint - The first value in the list is used to define the initial value of xval
    constr_func = ['Area']
    constr_values = [50]
    # Plot mesh  
    mesh_deform = True 
    factor = 10000
    # Save plots
    save = False
    #
    timing = False

    # interval = np.arange(freq_rsp[0], freq_rsp[1] + 1, freq_rsp[2])
    # func_vector = np.empty((len(interval), 5), dtype=complex)
    # origin = fc_opt.freq_resp(freq_rsp, const_func, constr_func, constr_values, force_matrix, restri_matrix, coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha_plot, beta_plot, eta_plot, x_min, p_par, q_par, func_name, None,None)
    # for i, modes in enumerate([1,2,3]):
    #     func_vector[:, i] = fc_opt.freq_resp(freq_rsp, const_func, constr_func, constr_values, force_matrix, restri_matrix, coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha_plot, beta_plot, eta_plot, x_min, p_par, q_par, func_name, modes, save=None)
        
    #     fig, ax = plt.subplots()
    #     ax.plot(interval, origin.real, label='original')
    #     ax.plot(interval, func_vector[:, i].real, label=str(modes) + 'modes')
    #     ax.set_xlabel('frequency [Hz]', fontsize=16)
    #     ax.set_ylabel(func_name.lower(), fontsize=16)
    #     ax.set_yscale('log')
    #     plt.legend()
    #     plt.savefig(str(modes) + ".eps")

    # fig, ax = plt.subplots()
    # ax.plot(interval, origin.real, label='original')
    # for i in range(3):
    #     ax.plot(interval, func_vector[:, i].real, label=str(modes) + 'modes')
    # ax.set_xlabel('frequency [Hz]', fontsize=16)
    # ax.set_ylabel(func_name.lower(), fontsize=16)
    # ax.set_yscale('log')
    # plt.legend()
    # plt.savefig('todospequenos' + ".eps")   
    #plt.show()


    beam.main(nelx, nely, lx, ly, func_name, force_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, freq_rsp, dens_filter, each_iter, mesh_deform=mesh_deform, factor=factor, max_iter=max_iter, save=save, timing=timing)