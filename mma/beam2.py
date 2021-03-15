from __future__ import division
from time import time
import numpy as np
import matplotlib.pyplot as plt
import pyqtgraph as pg
import logging
import sys
import os
sys.path.append(os.getcwd())
import solver_fem_2d.functions_2d as fc
import solver_fem_2d.plots_2d as plt_fem
import functions as opt
import plots as plt_opt
# Import MMA functions
from mma_opt import mmasub, subsolv, kktcheck

def main(nelx, nely, lx, ly, func_name, load_matrix, restri_matrix=None, freq1=180, constr_func=['Area'], constr_values=[50], n1=1, multiobjective=(None, 0), const_func=100, fac_ratio=2.1, modes=None, rho=7860, E=210e9, v=0.3, x_min=0.001, alpha_par=0, beta_par=5e-6, eta_par=0, alpha_plot=0, beta_plot=1e-8, eta_plot=0, p_par=3, q_par=1, freq_rsp=[], dens_filter=True, each_iter=True, max_iter=100, mesh_deform=False, factor=1000, save=False, timing=False):
    """ 
    Args:
        nelx (:obj:`int`): Number of elements on the x-axis.
        nely (:obj:`int`): Number of elements on the y-axis.
        lx (:obj:`int`): x-axis length.
        ly (:obj:`int`): x-axis length.
        func_name (:obj:`str`): Objective function used.
            It can be: 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'.
            If the multiobjective function is being calculated, weight n1 is assigned.
        load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value.
        restri_matrix (:obj:`numpy.array`, optional): The columns are respectively node, x direction, y direction. Defaults to None. 
        freq1 (:obj:`int`, optional): Optimized frequency. Defaults to 180.
        constr_func (:obj:`list`, optional): Restriction functions applied. Defaults to 'Area'.
            It can be: 'Area' or 'R Ratio'.
            The first function in the list will be used to define the initial value of xval.
        constr_values (:obj:`list`, optional): Values of restriction functions applied. Defaults to 50.
            Value in position i relates to the function in position i of the list constr_func.
            It can be a maximum of 4 values.
        n1 (:obj:`float`, optional): Weight n1 used in func_name. Defaults to 1.
            If n1 < 0: Maximize objective function
            If n1 > 0: Minimize objective function
        multiobjective (:obj:`tuple`, optional): Second function and frequency in the multiobjective. Defaults to (None, 0). 
            First value is the second function of the multiobjective function. The assigned weight is (1 - n1).
            Second value is the frequency that func_name2 is being optimized.            
        const_func (:obj:`float`, optional): Defaults to 100. 
        fac_ratio(:obj:`float`, optional): Factor applied in the radius to get elements in the vicinity of each element. Defaults to 2.1.
        modes (:obj:`int`, optional): If not None is used the Mode Superposition Method to calculate the displacement. Defaults to None.
        rho (:obj:`float`, optional): Density. Defaults to 7860. 
        E (:obj:`float`, optional): Elastic modulus. Defaults to 210e9.
        v (:obj:`float`, optional): Poisson's ratio. Defaults to 0.3. 
        x_min (:obj:`float`, optional): Minimum relative densities. Defaults to 0.001.
        alpha_par (:obj:`float`, optional): Damping coefficient proportional to mass. Defaults to 0.
        beta_par (:obj:`float`, optional): Damping coefficient proportional to stiffness. Defaults to 5e-6. 
        eta_par (:obj:`float`, optional): Damping coefficient. Defaults to 0.
        alpha_plot (:obj:`float`, optional): Damping coefficient proportional to mass to generate the graph. Defaults to 0. 
        beta_plot (:obj:`float`, optional): Damping coefficient proportional to stiffness to generate the graph. Defaults to 1e-8. 
        eta_plot (:obj:`float`, optional): Damping coefficient to generate the graph. Defaults to 0.
        p_par (:obj:`int`, optional): Penalization power to stiffness. Defaults to 3.
        q_par (:obj:`int`, optional): Penalization power to mass. Defaults to 1. 
        freq_rsp (:obj:`list`, optional): If len is 3, a frequency response graph of the original and optimized structure is generated. Defaults to [].
            First value is the minimum frequency of the graph.
            Second value is the maximum frequency of the graph.
            Third value is the step between each calculation of the objective function. 
        dens_filter (:obj:`bool`, optional): If True use density filter and False use sensitivity filter. Defaults to True.
        each_iter (:obj:`bool`, optional): If True plots the convergence graph for each iteration of the optimization. Defaults to True. 
        max_iter (:obj:`int`, optional): Number of iterations. Defaults to 100. 
        mesh_deform (:obj:`bool`, optional): If True plots the mesh deformation of the dynamic function. Defaults to False.
        factor (:obj:`float`, optional): Factor to deform the mesh. Defaults to 1000.
        save (:obj:`bool`, optional): if True save the optimization and frequency response graphs as PNG. Defaults to False.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.

    """
    t0 = time()
    # FEM settings
    coord, connect, ind_rows, ind_cols = fc.regularmeshQ4(lx, ly, nelx, nely)
    func_name2, freq2 = multiobjective
    omega1_par = 2 * np.pi * freq1
    omega2_par = 2 * np.pi * freq2   
    ngl = 2 * ((nelx + 1) * (nely + 1))
    # Beam initial settings
    m = len(constr_func)
    n = nelx * nely
    fval, dfdx = np.zeros((m, 1)), np.zeros((m, n))
    epsimin = 0.0000001
    eeen = np.ones((n,1))
    eeem = np.ones((m,1))
    zeron = np.zeros((n,1))
    zerom = np.zeros((m,1))
    initial_xval = opt.set_initxval(constr_func, constr_values)
    xval = initial_xval * eeen
    xval_original = xval
    xold1 = xval.copy()
    xold2 = xval.copy()
    xmin = 0.00 * eeen.copy()
    xmax = 1 * eeen
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
    # Get free indexes  
    free_ind = None
    if restri_matrix is not None:
        restricted_ind = fc.get_dofs(restri_matrix)
        free_ind = fc.remove_dofs(nelx, nely, restricted_ind)
    # Calculate load
    load_vector = fc.get_load_vector(nelx, nely, load_matrix)
    # Calculate neighbors and distance
    radius = fac_ratio * lx/nelx
    neighbors, H = opt.get_neighbors_radius(nelx, nely, coord, connect, radius)	
    # Calculate area and derivative
    area = opt.calc_A(coord, connect[:, 1:] - 1)
    # Calculate function values and gradients of the objective and constraints functions
    if outeriter == 0:   
        if dens_filter:
            xnew = opt.calc_xnew(H, neighbors, xval)
        else:
            xnew = xval
        stif_matrix, mass_matrix, dyna_stif, t_assembly = opt.solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_par, beta_par, eta_par, omega1_par, xnew, x_min, p_par, q_par)
        if modes is not None:
            disp_vector, t_superp = opt.mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega1_par, alpha_par, beta_par, eta_par, free_ind)
        else: 
            disp_vector, t_harmonic = opt.harmonic_problem(ngl, dyna_stif, load_vector, free_ind)
        # Area function
        fval, dfdx = opt.apply_constr(fval, dfdx, constr_func, constr_values, nelx, nely, lx, ly, coord, connect, E, v, rho, alpha_par, beta_par, p_par, q_par, x_min, area, xnew, disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, const_func, free_ind)
        # Objective function      
        if func_name == "Compliance":
            f0val = opt.func_compliance(disp_vector, load_vector)
            lam_par = opt.lambda_parameter(disp_vector, load_vector, f0val)
            df0dx = opt.derivative_compliance(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par)
        #
        elif func_name == "Elastic Potential Energy":
            f0val, fvirg = opt.elastic_potential_energy(disp_vector, stif_matrix, const_func)
            lam_par = opt.lambda_parameter_ep(disp_vector, stif_matrix, dyna_stif, free_ind)
            df0dx = opt.derivative_ep(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par, fvirg)
        #
        elif func_name == "Input Power":
            f0val, fvirg = opt.func_input_power(disp_vector, load_vector, omega1_par, const_func)
            df0dx = opt.derivative_input_power(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, fvirg)
        #               
        elif func_name == "Kinetic Energy":
            f0val, fvirg = opt.kinetic_energy(disp_vector, mass_matrix, omega1_par, const_func)
            lam_par = opt.lambda_parameter_ek(disp_vector, mass_matrix, dyna_stif, omega1_par, free_ind)
            df0dx = opt.derivative_ek(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par, fvirg)
        #   
        elif func_name == 'R Ratio':
            f0val, fvirg, elastic_p, kinetic_e = opt.R_ratio(disp_vector, stif_matrix, mass_matrix, omega1_par, const_func)
            lam_par = opt.lambda_parameter_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, fvirg, kinetic_e, free_ind)
            df0dx = opt.derivative_R(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par, fvirg, elastic_p, kinetic_e)
        # Normalization
        f0_scale = f0val
        f0val = n1 * 100 * f0val/f0_scale
        df0dx = n1 * 100 * df0dx/f0_scale
        # Multiobjective
        if (func_name2 is not None) and (n1 != 1):
            stif_matrix, mass_matrix, dyna_stif, t_assembly = opt.solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_par, beta_par, eta_par, omega2_par, xnew, x_min, p_par, q_par)
            if modes is not None:
                disp_vector2 = opt.mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega2_par, alpha_par, beta_par, eta_par, free_ind)
            else: 
                disp_vector2, t_harmonic = opt.harmonic_problem(ngl, dyna_stif, load_vector, free_ind)
            
            if func_name2 == "Compliance":
                f0val2 = opt.func_compliance(disp_vector2, load_vector)
                lam_par = opt.lambda_parameter(disp_vector2, load_vector, f0val2)
                df0dx2 = opt.derivative_compliance(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par)
            #
            elif func_name2 == "Elastic Potential Energy":
                f0val2, fvirg = opt.elastic_potential_energy(disp_vector2, stif_matrix, const_func)
                lam_par = opt.lambda_parameter_ep(disp_vector2, stif_matrix, dyna_stif, free_ind)
                df0dx2 = opt.derivative_ep(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par, fvirg)
            #
            elif func_name2 == "Input Power":
                f0val2, fvirg = opt.func_input_power(disp_vector2, load_vector, omega2_par, const_func)
                df0dx2 = opt.derivative_input_power(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, fvirg)
            #               
            elif func_name2 == "Kinetic Energy":
                f0val2, fvirg = opt.kinetic_energy(disp_vector2, mass_matrix, omega2_par, const_func)
                lam_par = opt.lambda_parameter_ek(disp_vector2, mass_matrix, dyna_stif, omega2_par, free_ind)
                df0dx2 = opt.derivative_ek(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par, fvirg)
            #   
            elif func_name2 == 'R Ratio':
                f0val2, fvirg, elastic_p, kinetic_e = opt.R_ratio(disp_vector2, stif_matrix, mass_matrix, omega2_par, const_func)
                lam_par = opt.lambda_parameter_R(disp_vector2, dyna_stif, stif_matrix, mass_matrix, omega2_par, fvirg, kinetic_e, free_ind)
                df0dx2 = opt.derivative_R(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par, fvirg, elastic_p, kinetic_e)
            # Normalization
            f0_scale_n2 = f0val2
            f0val2 = (1 - abs(n1)) * 100 * f0val2/f0_scale_n2
            df0dx2 = (1 - abs(n1)) * 100 * df0dx2/f0_scale_n2
            # Sum of functions and derivatives
            f0val = f0val + f0val2
            df0dx = df0dx + df0dx2
        # Apply density filter
        if dens_filter:
            dfdx = opt.dens_dconstr(dfdx, constr_func, H, neighbors, radius)
            df0dx = opt.density_filter(df0dx, H, neighbors) 
        else:
            dfdx = opt.sens_dconst(dfdx, constr_func, H, neighbors, xval, radius)
            df0dx = opt.sensitivity_filter(df0dx, H, neighbors, xval, radius)
            
        innerit = 0
        # Log
        set_logger(logger, outeriter, innerit, f0val, fval, xval)
    # List to plot convergence
    list_iter = []
    list_f0val = []
    list_fvals = []
    for i in range(m):
        list_fvals.append([])
    list_iter, list_f0val, list_fvals = opt.update_lists(outit, fval, f0val, list_iter, list_fvals, list_f0val, constr_func, constr_values)
    # Construct a QApplication
    app = pg.mkQApp()
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    if each_iter:
        win, p2, image = plt_opt.window_each_iter(constr_func, list_iter, list_f0val, list_fvals, xval, nelx, nely, func_name)
    else:
        gv, image = plt_opt.simple_window(xval, nelx, nely)
    #
    kktnorm = kkttol+10
    while (kktnorm > kkttol) and (outit < max_iter):
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
        stif_matrix, mass_matrix, dyna_stif, t_assembly = opt.solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_par, beta_par, eta_par, omega1_par, xnew, x_min, p_par, q_par)
        if modes is not None:
            disp_vector, t_superp = opt.mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega1_par, alpha_par, beta_par, eta_par, free_ind)
        else: 
            disp_vector, t_harmonic = opt.harmonic_problem(ngl, dyna_stif, load_vector, free_ind)
        # Area function
        fval, dfdx = opt.apply_constr(fval, dfdx, constr_func, constr_values, nelx, nely, lx, ly, coord, connect, E, v, rho, alpha_par, beta_par, p_par, q_par, x_min, area, xnew, disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, const_func, free_ind)
        # Objective function 
        if func_name == "Compliance":
            f0val = opt.func_compliance(disp_vector, load_vector)
            lam_par = opt.lambda_parameter(disp_vector, load_vector, f0val)
            df0dx = opt.derivative_compliance(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par)
        #
        elif func_name == "Elastic Potential Energy":
            f0val, fvirg = opt.elastic_potential_energy(disp_vector, stif_matrix, const_func)
            lam_par = opt.lambda_parameter_ep(disp_vector, stif_matrix, dyna_stif, free_ind)
            df0dx = opt.derivative_ep(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par, fvirg)
        #
        elif func_name == "Input Power":
            f0val, fvirg = opt.func_input_power(disp_vector, load_vector, omega1_par, const_func)
            df0dx = opt.derivative_input_power(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, fvirg)
        #               
        elif func_name == "Kinetic Energy":
            f0val, fvirg = opt.kinetic_energy(disp_vector, mass_matrix, omega1_par, const_func)
            lam_par = opt.lambda_parameter_ek(disp_vector, mass_matrix, dyna_stif, omega1_par, free_ind)
            df0dx = opt.derivative_ek(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par, fvirg)
        #   
        elif func_name == 'R Ratio':
            f0val, fvirg, elastic_p, kinetic_e = opt.R_ratio(disp_vector, stif_matrix, mass_matrix, omega1_par, const_func)
            lam_par = opt.lambda_parameter_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, fvirg, kinetic_e, free_ind)
            df0dx = opt.derivative_R(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min, xnew, disp_vector, lam_par, fvirg, elastic_p, kinetic_e)
        # Normalization
        f0val = n1 * 100 * f0val/f0_scale
        df0dx = n1 * 100 * df0dx/f0_scale
        # Multiobjective
        if (func_name2 is not None) and (n1 != 1):
            stif_matrix, mass_matrix, dyna_stif, t_assembly = opt.solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_par, beta_par, eta_par, omega2_par, xnew, x_min, p_par, q_par)
            if modes is not None:
                disp_vector2 = opt.mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega2_par, alpha_par, beta_par, eta_par, free_ind)
            else: 
                disp_vector2, t_harmonic = opt.harmonic_problem(ngl, dyna_stif, load_vector, free_ind)
            if func_name2 == "Compliance":
                f0val2 = opt.func_compliance(disp_vector2, load_vector)
                lam_par = opt.lambda_parameter(disp_vector2, load_vector, f0val2)
                df0dx2 = opt.derivative_compliance(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par)
            #
            elif func_name2 == "Elastic Potential Energy":
                f0val2, fvirg = opt.elastic_potential_energy(disp_vector2, stif_matrix, const_func)
                lam_par = opt.lambda_parameter_ep(disp_vector2, stif_matrix, dyna_stif, free_ind)
                df0dx2 = opt.derivative_ep(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par, fvirg)
            #
            elif func_name2 == "Input Power":
                f0val2, fvirg = opt.func_input_power(disp_vector2, load_vector, omega2_par, const_func)
                df0dx2 = opt.derivative_input_power(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, fvirg)
            #               
            elif func_name2 == "Kinetic Energy":
                f0val2, fvirg = opt.kinetic_energy(disp_vector2, mass_matrix, omega2_par, const_func)
                lam_par = opt.lambda_parameter_ek(disp_vector2, mass_matrix, dyna_stif, omega2_par, free_ind)
                df0dx2 = opt.derivative_ek(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par, fvirg)
            #   
            elif func_name2 == 'R Ratio':
                f0val2, fvirg, elastic_p, kinetic_e = opt.R_ratio(disp_vector2, stif_matrix, mass_matrix, omega2_par, const_func)
                lam_par = opt.lambda_parameter_R(disp_vector2, dyna_stif, stif_matrix, mass_matrix, omega2_par, fvirg, kinetic_e, free_ind)
                df0dx2 = opt.derivative_R(coord, connect, E, v, rho, alpha_par, beta_par, omega2_par, p_par, q_par, x_min, xnew, disp_vector2, lam_par, fvirg, elastic_p, kinetic_e)
            # Normalization
            f0val2 = (1 - abs(n1)) * 100 * f0val2/f0_scale_n2
            df0dx2 = (1 - abs(n1)) * 100 * df0dx2/f0_scale_n2
            # Sum of functions and derivatives
            f0val = f0val + f0val2
            df0dx = df0dx + df0dx2
        # Apply density filter
        if dens_filter:
            dfdx = opt.dens_dconstr(dfdx, constr_func, H, neighbors, radius)
            df0dx = opt.density_filter(df0dx, H, neighbors)
        else:
            dfdx = opt.sens_dconst(dfdx, constr_func, H, neighbors, xval, radius)
            df0dx = opt.sensitivity_filter(df0dx, H, neighbors, xval, radius)
        # The residual vector of the KKT conditions is calculated
        residu,kktnorm,residumax = \
            kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
        # Plot xval and objective function
        image.setImage(-xval.real.reshape((nely, nelx)))
        list_iter, list_f0val, list_fvals = opt.update_lists(outit, fval, f0val, list_iter, list_fvals, list_f0val, constr_func, constr_values)
        if each_iter:
            plt_opt.convergence(constr_func, p2, list_iter, list_f0val, list_fvals)
        pg.QtGui.QApplication.processEvents()
        # Log
        set_logger(logger, outeriter, innerit, f0val, fval, xval)
    # Final log
    logger.info(" Finished")     
    # Plot convergence
    if not each_iter:
        win, p2 = plt_opt.win_convergence(constr_func, list_iter, list_f0val, list_fvals, func_name)
    if save:
        plt_opt.save_fig(image, p2)
    tf = time()
    if timing:
        if modes is not None:
            print("Time to solve the Mode Superposition Method: " + t_superp + '[s]')
        else:
            print("Time to solve the harmonic analysis problem: " + t_harmonic + '[s]')
        print("Time to assembly global matrices: " + t_assembly + '[s]')
        print("Time to process: " + str(round((tf - t0), 6)) + '[s]')
    
    if len(freq_rsp) == 3:
        freq_range = freq_rsp[:2]
        delta = freq_rsp[2]
        f_original = opt.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_plot, beta_plot, eta_plot, xval_original, x_min, p_par, q_par, freq_range, delta, func_name, const_func, modes, load_vector, unrestricted_ind=free_ind)
        f_optimized = opt.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_plot, beta_plot, eta_plot, xval, x_min, p_par, q_par, freq_range, delta, func_name, const_func, modes, load_vector, unrestricted_ind=free_ind)
        ax = plt_opt.compare_freqresponse(freq_range, delta, f_optimized.real, f_original.real, func_name, save=save)
    if mesh_deform:
        disp_vector = fc.change_U_shape(disp_vector.real)
        coord_U = fc.apply_U(disp_vector, coord, factor)
        collection = plt_fem.build_collection(coord_U, connect[:, 1:])
        ax1 = plt_fem.plot_collection(lx, ly, coord_U, collection, load_matrix, restri_matrix, save=save)
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

def set_logger(logger, outeriter, innerit, f0val, fval, xval):
    outvector1 = np.array([outeriter, innerit, f0val, fval.flatten()])
    outvector2 = xval.flatten()
    # Log
    logger.info("outvector1 = {}".format(outvector1))
    logger.info("outvector2 = {}\n".format(outvector2))