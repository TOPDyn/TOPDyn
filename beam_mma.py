from __future__ import division
from time import time
import os
import logging
import numpy as np
import pyqtgraph as pg
import matplotlib.pyplot as plt

import functions_2d as fc
import plots_2d as plt_fem
import functions_opt as opt
import functions_obj as obj
import functions_deriv as df
import plots_opt as plt_opt
from mesh_process_2d import import_mesh
from mma_opt import mmasub, kktcheck

def main(mesh_file, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix=None, freq1=180, constr_func=["area"], constr_values=[50], n1=1, multiobjective=(None, 0), const_func=100, fac_ratio=2.1, modes=None, rho=7860, E=210e9, v=0.3, x_min_m=0.001, x_min_k=0.001, alpha_par=0, beta_par=5e-6, eta_par=0, alpha_plot=0, beta_plot=1e-8, eta_plot=0, p_par=3, q_par=1, passive_coord=None, freq_rsp=[], dens_filter=True, each_iter=True, max_iter=100, mesh_deform=False, factor=1000, save=False, timing=False):
    """ 
    Args:
        nelx (:obj:`int`): Number of elements on the x-axis.
        nely (:obj:`int`): Number of elements on the y-axis.
        lx (:obj:`float`): x-axis length.
        ly (:obj:`float`): x-axis length.
        func_name (:obj:`str`): Objective function used.
            
            * It can be "compliance", "input_power", "elastic_potential_energy", "kinetic_energy", "r_ratio", "local_ep", "local_ki" or "local_r".
            * The objective functions "local_ep", "local_ki" and "local_r" are designed to passive elements.
            * If the multiobjective function is being calculated, weight n1 is assigned.
        
        load_matrix (:obj:`numpy.array`): List of dictionaries. The dictionary can be:   
        
            * {"coord":value_coordinate, "axis":column_to_compare, "eps":error_margin, "x_direc":force_applied_x, "y_direc":force_applied_y, "force":force_value}
            * {"x_coord":x_coordinate, "y_coord":y_coordinate, "apply_x":force_applied_x, "apply_y":force_applied_y, "force":force_value}
                
        restri_matrix (:obj:`numpy.array`, optional): List of dictionaries. Defaults to None. The dictionary can be:
                
                * {"x_coord":x_coordinate, "y_coord":y_coordinate, "constrain_disp_x":constrain_disp_x, "constrain_disp_y":constrain_disp_y}
                * {"coord":value_coordinate, "axis":column_to_compare, "eps":error_margin, "constrain_disp_x":constrain_disp_x, "constrain_disp_y":constrain_disp_y} 
        
        freq1 (:obj:`int`, optional): Optimized frequency. Defaults to 180.
        constr_func (:obj:`list`, optional): Constraint functions applied. Defaults to "area".
            
            * It can be: "area", "r_ratio", "compliance", "local_ep", "local_ki" or "local_r".
            * The first function in the list will be used to define the initial value of xval.
            * If the same function is passed 2x,the box constraint is used. Negative values indicate the lower constraint.
            * Example: 
                * constr_func = ["area", "area"]
                * constr_values = [50, -20]
        
        constr_values (:obj:`list`, optional): Values of constraint functions applied. Defaults to 50.
                
                * Value in position i relates to the function in position i of the list constr_func.
                * constr_values[i] < 0 = lower constraint
                * constr_values[i] > 0 = upper constraint
                * If "compliance", "local_ep", "local_ki" or "local_r" is passed a tuple with constraint value and frequency respectively.
                * Example: 
                    * constr_func   = ["area", "area", "compliance", "r_ratio"]
                    * constr_values = [50, -20, (50, 1000), 10]
        
        n1 (:obj:`float`, optional): Weight n1 used in func_name. Defaults to 1.
            
            * If n1 < 0: Maximize objective function
            * If n1 > 0: Minimize objective function
        
        multiobjective (:obj:`tuple`, optional): Second function and frequency in the multiobjective. Defaults to (None, 0). 
            
            * First value is the second function of the multiobjective function. The assigned weight is (1 - n1).
            * Second value is the frequency that func_name2 is being optimized.            
        
        const_func (:obj:`float`, optional): Defaults to 100. 
        fac_ratio(:obj:`float`, optional): Factor applied in the radius to get elements in the vicinity of each element. Defaults to 2.1.
        modes (:obj:`int`, optional): If not None is used the Mode Superposition Method to calculate the displacement. Defaults to None.
        rho (:obj:`float`, optional): Density. Defaults to 7860. 
        E (:obj:`float`, optional): Elastic modulus. Defaults to 210e9.
        v (:obj:`float`, optional): Poisson's ratio. Defaults to 0.3. 
        x_min_m (:obj:`float`, optional): Minimum relative densities to mass. Defaults to 0.001.
        x_min_k (:obj:`float`, optional): Minimum relative densities to stiffness. Defaults to 0.001.
        alpha_par (:obj:`float`, optional): Damping coefficient proportional to mass. Defaults to 0.
        beta_par (:obj:`float`, optional): Damping coefficient proportional to stiffness. Defaults to 5e-6. 
        eta_par (:obj:`float`, optional): Damping coefficient. Defaults to 0.
        alpha_plot (:obj:`float`, optional): Damping coefficient proportional to mass to generate the graph. Defaults to 0. 
        beta_plot (:obj:`float`, optional): Damping coefficient proportional to stiffness to generate the graph. Defaults to 1e-8. 
        eta_plot (:obj:`float`, optional): Damping coefficient to generate the graph. Defaults to 0.
        p_par (:obj:`int`, optional): Penalization power to stiffness. Defaults to 3.
        q_par (:obj:`int`, optional): Penalization power to mass. Defaults to 1. 
        passive_coord (:obj:`tuple`): Region that the shape will not be changed. Defaults to None. 
                
                * Example: ((0.5, 1), (0.3, 0.6)) = ((x_initial, x_final), (y_initial, y_final))
        
        freq_rsp (:obj:`list`, optional): If len is 3, a frequency response graph of the original and optimized structure is generated. Defaults to [].
            
            * First value is the minimum frequency of the graph.
            * Second value is the maximum frequency of the graph.
            * Third value is the step between each calculation of the objective function. 
        
        dens_filter (:obj:`bool`, optional): If True use density filter and False use sensitivity filter. Defaults to True.
        each_iter (:obj:`bool`, optional): If True plots the convergence graph for each iteration of the optimization. Defaults to True. 
        max_iter (:obj:`int`, optional): Number of iterations. Defaults to 100. 
        mesh_deform (:obj:`bool`, optional): If True plots the mesh deformation of the dynamic function. Defaults to False.
        factor (:obj:`float`, optional): Factor to deform the mesh. Defaults to 1000.
        save (:obj:`bool`, optional): if True save the optimization and frequency response graphs as PNG. Defaults to False.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.
    """
    # Check function names
    opt.check_func_name(func_name, multiobjective[0], constr_func)
    
    t0 = time()
    # FEM settings
    if mesh_file is not None:
        path = os.path.dirname(os.path.realpath(__file__)) 
        m_file = os.path.join(path, mesh_file)
        coord, connect = import_mesh(m_file)
        ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)
        nelx = len(coord[coord[:, 2] == coord[0, 2]]) - 1
        nely = len(coord[coord[:, 1] == coord[0, 1]]) - 1
        lx = max(coord[:, 1])
        ly = max(coord[:, 2])
    else:
        coord, connect = fc.regularmeshQ4(lx, ly, nelx, nely, timing=timing)
    ind_dofs = fc.get_ind_dofs(connect, 2)
    ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)

    # Load matrix
    load_matrix = fc.get_matrices(load_matrix, coord, True)

    # Get free indexes  
    free_ind = None
    if restri_matrix is not None:
        restri_matrix = fc.get_matrices(restri_matrix, coord, False)
        restricted_ind = fc.get_dofs(restri_matrix)
        free_ind = fc.remove_dofs(nelx, nely, restricted_ind)
    
    func_name2, freq2 = multiobjective
    multiobj_bool = (func_name2 is not None) and (n1 != 1)
    omega1_par = 2 * np.pi * freq1
    omega2_par = 2 * np.pi * freq2   
    ngl = 2 * ((nelx + 1) * (nely + 1))
    natural_freqs = None

    # Constrain settings
    aux = ["compliance", "local_ep", "local_ki", "local_r"]
    freq_constr_bool = any(x in aux for x in constr_func)
    if freq_constr_bool:
        constr_values, freq_constr, ind_freq_constr, ind_constr2 = opt.constr_with_freq(constr_func, constr_values)
        f_scale_constr = 100*np.ones(len(constr_values))
    else:
        freq_constr = None
        ind_constr2 = None
    
    constr_values = np.array(constr_values)

    # Calculate neighbors and distance
    radius = fac_ratio * lx/nelx
    neighbors, H, centroids = opt.get_neighbors_radius(nelx, nely, coord, connect, radius)

    # Beam initial settings
    m = len(constr_func)
    n = nelx * nely
    fval, dfdx = np.zeros((m, 1)), np.zeros((m, n))
    eeen = np.ones((n,1))
    eeem = np.ones((m,1))
    zerom = np.zeros((m,1))
    initial_xval = opt.set_initxval(constr_func, constr_values)
    xval = initial_xval * eeen
    xval_original = xval
    xold1 = xval.copy()
    xold2 = xval.copy()
    xmin = 0.00 * eeen.copy()
    xmax = 1 * eeen
    
    if passive_coord is not None:
        passive_el = opt.get_passive_el(passive_coord, centroids)
        xmin, xval = opt.set_passive_el(xmin, xval, passive_el)
        ind_passive = ind_dofs[passive_el, :]
    else:
        passive_el = None
        ind_passive = None

    low = xmin.copy()
    upp = xmax.copy()
    move = 1.0
    c = 1000*eeem 
    d = eeem.copy()
    a0 = 1
    a = zerom.copy()
    outeriter = 0 
    kkttol = 0	
    outit = 0
    # Logger
    path = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(path, "MMA_BEAM2.log")
    logger = setup_logger(file)

    # Calculate load
    load_vector = fc.get_load_vector(nelx, nely, load_matrix)

    # Calculate area and derivative
    area = opt.calc_A(coord, connect[:, 1:] - 1)

    if save:
        folder_name = 'data'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        os.makedirs(directory, exist_ok=True)

        xnew_dir = os.path.join(directory, 'xnew')
        os.makedirs(xnew_dir, exist_ok=True)

        save_fval = np.empty((max_iter + 1, len(constr_func)))
        save_f0val1 = np.empty(max_iter + 1)
        save_fvirg1 = np.empty(max_iter + 1)
    
    # Calculate function values and gradients of the objective and constraints functions
    if outeriter == 0:   
        if dens_filter:
            xnew = opt.calc_xnew(H, neighbors, xval)
        else:
            xnew = xval
        data_k, data_m, _ = opt.solution2D(coord, connect, nelx, nely, E, v, rho, xnew, x_min_m, x_min_k, p_par, q_par)
        stif_matrix, mass_matrix, damp_matrix = opt.assembly_matrices(data_k, data_m, ind_rows, ind_cols, ngl, alpha_par, beta_par)
        dyna_stif = opt.assembly_dyna_stif(omega1_par, mass_matrix, damp_matrix, stif_matrix)
        disp_vector, natural_freqs, _ = opt.get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif, load_vector, free_ind, omega1_par, alpha_par, beta_par, eta_par, ngl)
        
        # Constraint function      
        fval, dfdx = opt.apply_constr(fval, dfdx, constr_func, constr_values, freq_constr, ind_constr2, lx, ly, ngl, coord, connect, E, v, rho, alpha_par, beta_par, eta_par, \
            p_par, q_par, x_min_m, x_min_k, area, xnew, modes, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, load_vector, omega1_par, const_func, free_ind, passive_el, ind_dofs)

        # Objective function      
        f0val, fvirg = obj.objective_funcs(func_name, disp_vector, stif_matrix, mass_matrix, load_vector, omega1_par, const_func, passive_el, ind_passive, coord, connect, E, v, rho)
        f0_scale = f0val
        # Derivative
        df0dx =  df.derivatives_objective(func_name, fvirg, disp_vector, coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min_m, x_min_k, xnew, \
                                        load_vector, mass_matrix, stif_matrix, dyna_stif, free_ind, ind_dofs, ngl, ind_passive, passive_el)

        if save:
            np.savetxt(os.path.join(xnew_dir, 'xnew'+'_'+str(outit)+'.txt'), xnew)
            save_fval[outit,:] = fval[:,0]              
            save_fvirg1[outit] = fvirg
            save_f0val1[outit] = f0val

        if freq_constr_bool:
            f_scale_constr[ind_freq_constr] = fval[ind_freq_constr, 0]
            fval = opt.normalize_constr(fval, f_scale_constr)
                
        fval[:, 0] -= constr_values

        # Multiobjective
        if multiobj_bool:
            dyna_stif2 = opt.assembly_dyna_stif(omega2_par, mass_matrix, damp_matrix, stif_matrix)
            disp_vector2, _, _ = opt.get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif2, load_vector, free_ind, omega2_par, alpha_par, beta_par, eta_par, ngl)
            
            # Second objective function
            f0val2, fvirg2 = obj.objective_funcs(func_name2, disp_vector2, stif_matrix, mass_matrix, load_vector, omega2_par, const_func, passive_el, ind_passive, coord, connect, E, v, rho)
            # Derivative
            df0dx2 = df.derivatives_objective(func_name2, fvirg2, disp_vector2, coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min_m, x_min_k, xnew, \
                                            load_vector, mass_matrix, stif_matrix, dyna_stif2, free_ind, ind_dofs, ngl, ind_passive, passive_el)

            if save:
                save_f0val2 = np.empty(max_iter + 1)
                save_fvirg2 = np.empty(max_iter + 1)
                save_f0val2[outit] = f0val2
                save_fvirg2[outit] = fvirg2
            
            # Filter
            if dens_filter:
                dfdx, df0dx, df0dx2 = opt.density_filter(H, neighbors, dfdx, df0dx, df0dx2)
            else:
                dfdx, df0dx, df0dx2 = opt.sensitivity_filter(H, neighbors, x_min_k, xnew, dfdx, df0dx, df0dx2)
            # Normalize multiobjective
            f0_scale_n2  = f0val2
            f0val2, df0dx2 = opt.normalize(1 - abs(n1), f0_scale_n2, f0val2, df0dx2)
            f0val, df0dx = opt.normalize(n1, f0_scale, f0val, df0dx)
            # Sum of functions and derivatives
            f0val = f0val + f0val2
            df0dx = df0dx + df0dx2
        else:
            # Filter
            if dens_filter:
                dfdx, df0dx, _ = opt.density_filter(H, neighbors, dfdx, df0dx)
            else:
                dfdx, df0dx, _ = opt.sensitivity_filter(H, neighbors, x_min_k, xnew, dfdx, df0dx)
            # Normalize objective f
            f0val, df0dx = opt.normalize(n1, f0_scale, f0val, df0dx)   
        innerit = 0
        # Log
        set_logger(logger, outeriter, innerit, f0val, fval, xval, natural_freqs)
    
    # List to plot convergence
    list_iter  = np.empty(max_iter + 1)
    list_f0val = np.empty(max_iter + 1)
    list_fvals = np.empty((max_iter + 1, len(constr_func)))
    list_iter, list_f0val, list_fvals = opt.update_lists(outit, fval, f0val, list_iter, list_fvals, list_f0val, constr_values)

    # Construct a QApplication
    app = pg.mkQApp()
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    labels_constr = plt_opt.legend_constr(constr_func)
    if each_iter:
        win, p2, curves_funcs, grid = plt_opt.window_each_iter(constr_func, func_name, labels_constr)
        plt_opt.set_conv_data(outeriter, curves_funcs, list_iter, list_f0val, list_fvals, constr_func)
    else:
        gv, grid = plt_opt.simple_window()
    x_plot, y_plot = plt_opt.set_coord_grid(lx, ly, nelx, nely)
    plt_opt.set_grid_data(grid, xnew, x_plot, y_plot, nelx, nely)
    pg.QtGui.QApplication.processEvents()
    
    kktnorm = kkttol+10
    chmax = 1
    chtol = 1e-4
    kconv = 0
    while (kktnorm > kkttol) and (outit < max_iter) and (kconv < 5):
        outit += 1
        outeriter += 1
        # The MMA subproblem is solved at the point xval:
        xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp = \
            mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,move)
        # Some vectors are updated:
        xold2 = xold1.copy()
        xold1 = xval.copy()
        xval = xmma.copy()

        # Re-calculate function values and gradients of the objective and constraints functions
        if dens_filter:
            xnew = opt.calc_xnew(H, neighbors, xval)
        else:
            xnew = xval
        data_k, data_m, t_assembly = opt.solution2D(coord, connect, nelx, nely, E, v, rho, xnew, x_min_m, x_min_k, p_par, q_par)
        stif_matrix, mass_matrix, damp_matrix = opt.assembly_matrices(data_k, data_m, ind_rows, ind_cols, ngl, alpha_par, beta_par)
        dyna_stif = opt.assembly_dyna_stif(omega1_par, mass_matrix, damp_matrix, stif_matrix)
        disp_vector, natural_freqs, t_U = opt.get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif, load_vector, free_ind, omega1_par, alpha_par, beta_par, eta_par, ngl)

        # Constraint function              
        fval, dfdx = opt.apply_constr(fval, dfdx, constr_func, constr_values, freq_constr, ind_constr2, lx, ly, ngl, coord, connect, E, v, rho, alpha_par, beta_par, eta_par, \
            p_par, q_par, x_min_m, x_min_k, area, xnew, modes, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, load_vector, omega1_par, const_func, free_ind, passive_el, ind_dofs)       

        # Objective function 
        f0val, fvirg = obj.objective_funcs(func_name, disp_vector, stif_matrix, mass_matrix, load_vector, omega1_par, const_func, passive_el, ind_passive, coord, connect, E, v, rho)
        # Derivative
        df0dx = df.derivatives_objective(func_name, fvirg, disp_vector, coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min_m, x_min_k, xnew, \
                                        load_vector, mass_matrix, stif_matrix, dyna_stif, free_ind, ind_dofs, ngl, ind_passive, passive_el)

        if save:
            np.savetxt(os.path.join(xnew_dir, 'xnew'+'_'+str(outit)+'.txt'), xnew)
            save_fval[outit,:] = fval[:,0]             
            save_fvirg1[outit] = fvirg
            save_f0val1[outit] = f0val

        if freq_constr_bool:
            fval = opt.normalize_constr(fval, f_scale_constr)
        
        fval[:, 0] -= constr_values
        
        # Multiobjective
        if multiobj_bool:
            dyna_stif2 = opt.assembly_dyna_stif(omega2_par, mass_matrix, damp_matrix, stif_matrix)
            disp_vector2, _, _ = opt.get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif2, load_vector, free_ind, omega2_par, alpha_par, beta_par, eta_par, ngl)

            # Second objective function
            f0val2, fvirg2 = obj.objective_funcs(func_name2, disp_vector2, stif_matrix, mass_matrix, load_vector, omega2_par, const_func, passive_el, ind_passive, coord, connect, E, v, rho)
            # Derivative
            df0dx2 = df.derivatives_objective(func_name2, fvirg2, disp_vector2, coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min_m, x_min_k, xnew, \
                                        load_vector, mass_matrix, stif_matrix, dyna_stif2, free_ind, ind_dofs, ngl, ind_passive, passive_el)
            
            if save:
                save_f0val2[outit] = f0val2
                save_fvirg2[outit] = fvirg2

            # Filter
            if dens_filter:
                dfdx, df0dx, df0dx2 = opt.density_filter(H, neighbors, dfdx, df0dx, df0dx2)
            else:
                dfdx, df0dx, df0dx2 = opt.sensitivity_filter(H, neighbors, x_min_k, xnew, dfdx, df0dx, df0dx2)

            # Normalize multiobjective
            f0val2, df0dx2 = opt.normalize(1 - abs(n1), f0_scale_n2, f0val2, df0dx2)
            f0val, df0dx = opt.normalize(n1, f0_scale, f0val, df0dx)
            # Sum of functions and derivatives
            f0val = f0val + f0val2
            df0dx = df0dx + df0dx2
        else:
            # Filter
            if dens_filter:
                dfdx, df0dx, _ = opt.density_filter(H, neighbors, dfdx, df0dx)
            else:
                dfdx, df0dx, _ = opt.sensitivity_filter(H, neighbors, x_min_k, xnew, dfdx, df0dx)
            # Normalize objective f
            f0val, df0dx = opt.normalize(n1, f0_scale, f0val, df0dx)   

        # The residual vector of the KKT conditions is calculated
        _, kktnorm, _ = \
            kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
        
        # conv. crit.
        if outeriter > 10: 
            chmax = max(abs(xold2 - xold1))/max(xold1)
            if chmax < chtol:
                kconv = kconv + 1
        else:
            chmax = 1

        # Plot xval and objective function   
        plt_opt.set_grid_data(grid, xnew, x_plot, y_plot, nelx, nely)
        list_iter, list_f0val, list_fvals = opt.update_lists(outit, fval, f0val, list_iter, list_fvals, list_f0val, constr_values)
        if each_iter:
            plt_opt.set_conv_data(outeriter, curves_funcs, list_iter, list_f0val, list_fvals, constr_func)
        pg.QtGui.QApplication.processEvents()
        # Log
        set_logger(logger, outeriter, innerit, f0val, fval, xval, natural_freqs)
    # Final log
    logger.info(" Finished")     
    # Plot convergence
    if not each_iter:
        win, p2 = plt_opt.win_convergence(constr_func, list_iter, list_f0val, list_fvals, func_name, labels_constr)

    tf = time()
    if timing:
        if modes is not None:
            print("Time to solve the Mode Superposition Method: " + t_U + '[s]')
        else:
            print("Time to solve the harmonic analysis problem: " + t_U + '[s]')
        print("Time to assembly global matrices: " + t_assembly + '[s]')
        print("Time to process: " + str(round((tf - t0), 6)) + '[s]')
    
    if len(freq_rsp) == 3:
        print("Calculating the frequency response of the objective function")
        print('initial conditions')
        f_original  = opt.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_plot, beta_plot, eta_plot, xval_original, x_min_m, x_min_k, p_par, q_par, freq_rsp, func_name, const_func, modes, load_vector, passive_el, ind_passive, aux_R=False, unrestricted_ind=free_ind)
        print('optimized conditions')
        f_optimized = opt.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_plot, beta_plot, eta_plot, xnew, x_min_m, x_min_k, p_par, q_par, freq_rsp, func_name, const_func, modes, load_vector, passive_el, ind_passive, aux_R=False, unrestricted_ind=free_ind)
        fig_freq, ax = plt_opt.compare_freqresponse(freq_rsp, f_optimized.real, f_original.real, func_name)
    
    if mesh_deform:
        disp_vector = fc.change_U_shape(disp_vector.real)
        coord_U = fc.apply_U(disp_vector, coord, factor)
        collection = plt_fem.build_collection(coord_U, connect[:, 1:])
        fig_mesh, ax1 = plt_fem.plot_collection(lx, ly, coord_U, collection, load_matrix, restri_matrix)

    if save:
        if multiobj_bool:
            ind_fval = 5
            aux_multi = 2
        else:
            ind_fval = 3
            aux_multi = 0
         
        data = np.empty((outit+1, 3+len(constr_func)+aux_multi))
        data[:, 0] = list_iter[:outit+1]

        data[:, 1] = save_f0val1[:outit+1]
        data[:, 2] = save_fvirg1[:outit+1]

        if multiobj_bool:
            data[:, 3] = save_f0val2[:outit+1]
            data[:, 4] = save_fvirg2[:outit+1]

        data[:, ind_fval:] = save_fval[:outit+1,:]

        header = opt.create_header(multiobj_bool, constr_func)
        np.savetxt(os.path.join(directory, 'functions.txt'), data, delimiter=",", header=header, comments='')          

        img_dir = os.path.join(directory, 'images')
        os.makedirs(img_dir, exist_ok=True)

        plt_opt.save_fig(grid, os.path.join(img_dir, 'opt.png'), pg_graph=True)
        plt_opt.save_fig(p2, os.path.join(img_dir, 'convergence.png'), pg_graph=True)

        if len(freq_rsp) == 3:
            np.savetxt(os.path.join(directory, 'frequency_rsp.txt'), np.column_stack((f_original, f_optimized)), delimiter=",", header="orig,opt", comments='')
            plt_opt.save_fig(fig_freq, os.path.join(img_dir, 'freqrsp.png'), pg_graph=False)
        
        if mesh_deform:
            plt_opt.save_fig(fig_mesh, os.path.join(img_dir, 'mesh.png'), pg_graph=False)

    print('Done!')
    plt.show()
    app.exec_()
 
# Setup logger
def setup_logger(logfile):
    # Create logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # Create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # Create file handler and set level to debug
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    # Add formatter to ch and fh
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # Add ch and fh to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
    # Set numpy print options
    np.set_printoptions(precision=4, formatter={'float': '{: 0.4f}'.format})
    # Open logfile and reset
    with open(logfile, 'w'): pass
    # Return logger
    return logger

def set_logger(logger, outeriter, innerit, f0val, fval, xval, natural_freqs):
    outvector1 = np.array([outeriter, innerit, f0val, fval.flatten()])
    outvector2 = xval.flatten()
    # Log
    logger.info("outvector1 = {}".format(outvector1))
    logger.info("outvector2 = {}\n".format(outvector2))
    if natural_freqs is not None:
        logger.info("Natural frequencies= {}\n".format(outvector2))