import cmath
import os
import numpy as np
from time import time
from scipy.linalg import eigh
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

import functions_deriv as df
import functions_2d as fc
import functions_obj as obj
import plots_opt as plt_opt
from mesh_process_2d import import_mesh

def exe_opt(mma, mesh_file, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix=None, freq1=180, constr_func=['Area'], constr_values=[50], n1=1, multiobjective=(None, 0), const_func=100, fac_ratio=2.1, modes=None, rho=7860, E=210e9, v=0.3, x_min_m=0.001, x_min_k=0.001, alpha_par=0, beta_par=5e-6, eta_par=0, alpha_plot=0, beta_plot=1e-8, eta_plot=0, p_par=3, q_par=1, passive_coord=None, freq_rsp=[], chtol=1e-4, dens_filter=True, each_iter=True, max_iter=100, mesh_deform=False, factor=1000, save=False, timing=False):
    if mma:
        import beam_mma as beam
        beam.main(mesh_file, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, passive_coord, freq_rsp, chtol, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)
    else:
        import beam_gcmma as beam
        beam.main(mesh_file, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, passive_coord, freq_rsp, chtol, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)

def solution2D(coord, connect, nelx, nely, E, v, rho, xval, x_min_m, x_min_k, p_par, q_par):
    """ Assembly matrices.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        ngl (:obj:`int`): Degrees of freedom.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.  
        rho (:obj:`float`): Density.  
        xval (:obj:`numpy.array`): Indicates where there is mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness.
        p_par (:obj:`int`): Penalization power to stiffness.
        q_par (:obj:`int`): Penalization power to mass. 

    Returns:
        A tuple with stiffnes, mass and time to assembly the matrices.
    """
    t01 = time() 
    data_k = np.zeros((nelx * nely, 64), dtype=complex)
    data_m = np.zeros((nelx * nely, 64), dtype=complex)
    
    for el in range(nelx * nely):
        Ke, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        data_k[el, :] = (x_min_k + (xval[el]**p_par)*(1-x_min_k))* Ke.flatten()
        if xval[el]>0.1:
            data_m[el, :] = (x_min_m + (xval[el]**q_par)*(1-x_min_m)) * Me.flatten()
        else:
            data_m[el, :] =  (x_min_m + (3.512e7*xval[el]**9 - 2.081e8*xval[el]**10)*(1-x_min_m) ) * Me.flatten() 
    
    tf1 = time()
    t_assembly = str(round((tf1 - t01), 6))
    return data_k.flatten(), data_m.flatten(), t_assembly
    
def assembly_matrices(data_k, data_m, ind_rows, ind_cols, ngl, alpha, beta):
    stif_matrix = csc_matrix((data_k, (ind_rows, ind_cols)), shape=(ngl, ngl))
    mass_matrix = csc_matrix((data_m, (ind_rows, ind_cols)), shape=(ngl, ngl))
    damp_matrix = alpha * mass_matrix + (beta) * stif_matrix
    return stif_matrix.real, mass_matrix.real, damp_matrix.real

def assembly_dyna_stif(omega_par, mass_matrix, damp_matrix, stif_matrix):
    return -(omega_par**2) * mass_matrix + 1j * omega_par * damp_matrix + stif_matrix

def harmonic_problem(ngl, dyna_stif, load_vector, free_ind=None):
    ''' Calculate displacement.

    Args:
        ngl (:obj:`int`): Degrees of freedom.
        dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix. 
        load_vector (:obj:`numpy.array`): Force.
        free_ind (:obj:`numpy.array`): Free dofs. 
        
    Returns:
        A tuple with displacement vector and time.
    '''
    t02 = time()
    if free_ind is not None:
        disp_vector = np.zeros((ngl), dtype=complex)
        disp_vector[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], load_vector[free_ind])
    else:
        disp_vector = spsolve(dyna_stif, load_vector)

    tf2 = time()
    t_harmonic = str(round((tf2 - t02),6))
    return disp_vector, t_harmonic

def modal_analysis(stif_matrix, mass_matrix, modes=20):
    """ Modal Analysis. Use eigh Scipy function.

    Args:
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        modes (:obj:`int`, optional): The number of eigenvalues and eigenvectors desired.
    
    Returns:
        A tuple with natural frequencies and modes shape.
    """
    eigen_values, eigen_vectors = eigh(stif_matrix.todense(), b=mass_matrix.todense(), subset_by_index=[0,modes])

    # mod_freqmax = 3 * freqmax
    # eigen_values, eigen_vectors = eigh(stif_matrix.todense(), b=mass_matrix.todense(), subset_by_value=[-np.inf, (2*np.pi*mod_freqmax)])

    positive_real = np.absolute(np.real(eigen_values))
    natural_frequencies = np.sqrt(positive_real)/(2*np.pi)
    modal_shape = np.real(eigen_vectors)

    index_order = np.argsort(natural_frequencies)
    natural_frequencies = natural_frequencies[index_order]
    modal_shape = modal_shape[:, index_order]
    return natural_frequencies, modal_shape

def mode_superposition(natural_frequencies, modal_shape, stif_matrix, load_vector, omega_par, alpha, beta, eta, free_ind):    
    """ Perform an harmonic analysis through superposition method and returns the response of
        all nodes due the external or internal equivalent load. It has been implemented two
        different damping models: Viscous Proportional and Hysteretic Proportional
        Entries for Viscous Proportional Model Damping: (alpha_v, beta_v)
        Entries for Hyteretic Proportional Model Damping: (alpha_h, beta_h)
    
    Args:
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        load_vector (:obj:`numpy.array`): Force.
        omega_par (:obj:`float`): 2 pi frequency
        alpha (:obj:`float`): Damping coefficient proportional to mass. 
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        eta (:obj:`float`): Damping coefficient. 
        free_ind (:obj:`numpy.array`): Free dofs. 

    Returns:
        A tuple with displacement and time to solve the problem.
    """
    t0 = time()
    alphaV, betaV, betaH = alpha, beta, eta

    #natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)

    F_aux = modal_shape.T @ load_vector[free_ind]
    omega_n = 2*np.pi*natural_frequencies
    F_kg = (omega_n**2)

    F_mg =  - (omega_par**2)
    F_cg = 1j*((betaH + betaV*omega_par)*(omega_n**2) + (0. + omega_par*alphaV)) 
    data = np.divide(1, (F_kg + F_mg + F_cg))
    diag = np.diag(data)
    #disp_vector = modal_shape @ (diag @ F_aux[:,i])
    rows = stif_matrix.shape[0]
    disp_vector = np.zeros((rows), dtype=complex)
    disp_vector[free_ind] = modal_shape @ (diag @ F_aux)

    tf = time()
    t_superp = str(round((tf - t0), 6))
    return disp_vector, t_superp

def get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif, load_vector, free_ind, omega_par, alpha_par, beta_par, eta_par, ngl):
    if modes is not None:
        natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)
        disp_vector, t_U = mode_superposition(natural_frequencies, modal_shape, stif_matrix, load_vector, omega_par, alpha_par, beta_par, eta_par, free_ind)
    else:  
        disp_vector, t_U = harmonic_problem(ngl, dyna_stif, load_vector, free_ind)
        natural_frequencies = None
    return disp_vector, natural_frequencies, t_U

def freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_m, x_min_k, p_par, q_par, freq_range, delta, func_name, const_func, modes, load_vector, passive_el, ind_passive, aux_R, **kwargs):
    """ Calculates the objective function for a range of frequencies.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        ind_rows (:obj:`numpy.array`): Node indexes to make the assembly.
        ind_cols (:obj:`numpy.array`): Node indexes to make the assembly.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        ngl (:obj:`int`): Degrees of freedom.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.  
        rho (:obj:`float`): Density.  
        alpha (:obj:`float`): Damping coefficient proportional to mass. 
        beta (:obj:`float`): Damping coefficient proportional to stiffness.  
        eta (:obj:`float`): Damping coefficient. 
        xval (:obj:`numpy.array`): Indicates where there is mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness. 
        p_par (:obj:`int`): Penalization power to stiffness. 
        q_par (:obj:`int`): Penalization power to mass.
        freq_range (:obj:`list`): Frequency range.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (:obj:`int`): Step between each calculation of the objective function. 
        func_name (:obj:`str`): Objective function used.
        const_func (:obj:`float`):
        modes (:obj:`int`): The number of eigenvalues and eigenvectors desired.
        load_vector (:obj:`numpy.array`): Force.

    Returns:
        Objective function values.
    """
    interval = np.arange(freq_range[0], freq_range[1] + 1, delta)
    func_vector = np.empty((len(interval)), dtype=complex)

    l = len(interval) + 5
    progress = 0
    printProgressBar(progress, l, prefix = 'Progress:', suffix = 'Complete', length = 50)

    free_ind = None
    if kwargs.get('unrestricted_ind') is not None:
        free_ind = kwargs.get('unrestricted_ind')

    data_k, data_m, _ = solution2D(coord, connect, nelx, nely, E, v, rho, xval, x_min_m, x_min_k, p_par, q_par)
    stif_matrix, mass_matrix, damp_matrix = assembly_matrices(data_k, data_m, ind_rows, ind_cols, ngl, alpha, beta)
    if modes is not None:
        natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)
    
    progress += 5
    printProgressBar(progress, l, prefix = 'Progress:', suffix = 'Complete', length = 50)

    for n in range(len(interval)):
        omega_par = 2 * np.pi * interval[n]
        dyna_stif = assembly_dyna_stif(omega_par, mass_matrix, damp_matrix, stif_matrix)
        if modes is not None:
            disp_vector, _ = mode_superposition(natural_frequencies, modal_shape, stif_matrix, load_vector, omega_par, alpha, beta, eta, free_ind)
        else: 
            disp_vector, _ = harmonic_problem(ngl, dyna_stif, load_vector, free_ind)

        _, func_vector[n] = obj.objective_funcs(func_name, disp_vector, stif_matrix, mass_matrix, load_vector, omega_par, const_func, passive_el, ind_passive, coord, connect, E, v, rho, aux_R=False)
        progress += 1
        printProgressBar(progress, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
    return func_vector 

def set_deriv(dfdx, df0dx, df0dx2):
    if df0dx2 is not None:
        cols = 2 + dfdx.shape[0]
    else:
        cols = 1 + dfdx.shape[0]
    
    all_deriv_f = np.empty((df0dx.shape[0], cols))
    for i in range(dfdx.shape[0]):
        all_deriv_f[:, i] = dfdx[i, :]
    
    if df0dx2 is not None:
        all_deriv_f[:, cols-1] = df0dx2[:, 0]
        all_deriv_f[:, cols-2] = df0dx[:, 0]
    else:
        all_deriv_f[:, cols-1] = df0dx[:, 0]
    return np.empty((df0dx.shape[0], cols)), all_deriv_f, cols

def out_deriv(all_deriv, dfdx, df0dx2):
    aux = dfdx.shape[0]
    if df0dx2 is None:
        return all_deriv[:, :aux].T, all_deriv[:, aux].reshape(-1, 1), None
    else:
        return all_deriv[:, :aux].T, all_deriv[:, aux].reshape(-1, 1), all_deriv[:, aux+1].reshape(-1, 1)

def get_neighbors_radius(nelx, nely, coord, connect, radius):
    """ Check neighboring elements that have the centroid within the predetermined radius.

    Args:
        nelx (:obj:`int`): Number of elements on the x axis.
        nely (:obj:`int`): Number of elements on the x axis
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        radius (:obj:`float`): Radius to get elements in the vicinity of each element.

    Returns:
        neighbors, H, centroids 
    """
    el_number = nelx * nely
     
    centroids = np.empty((el_number, 2))
    idx = connect[:, 1:] - 1   
    centroids[:, 0] = np.sum(coord[idx, 1], axis = 1)/4
    centroids[:, 1] = np.sum(coord[idx, 2], axis = 1)/4
     
    ind_rows = []
    ind_cols = []
    data = []
    cols = 0
    neighbors = []
    for el in range(el_number):
        distance = np.sqrt(np.sum((centroids[el] - centroids)**2, axis=1))
        mask = distance <= radius
        neighbor = mask.nonzero()[0] + 1
        neighbors.extend(neighbor - 1)
        
        hi      = radius - distance
        hi_max  = np.maximum(0, hi)
        data.extend(hi_max[mask])
        aux     = len(hi_max[mask])
        rows    = np.repeat(el, aux) #.tolist()
        columns = np.arange(0, aux)
        ind_rows.extend(rows) 
        ind_cols.extend(columns)
        
        if aux > cols:
            cols = aux
    H = csc_matrix((data, (ind_rows, ind_cols)), shape=(nelx*nely, cols)).toarray()
    neighbors = csc_matrix((neighbors, (ind_rows, ind_cols)), shape=(nelx*nely, cols), dtype='int').toarray()
    return neighbors, H, centroids

def density_filter(H1, neighbors, dfdx, df0dx, df0dx2=None):
    """ Apply the density filter to the derivative of the functions (constrain, objective and multiobjective).

    Args:
        H (:obj:`csc_matrix`): Radius subtracted from the distance between the element and the neighbors.
        neighbors (:obj:`csc_matrix`): Neighbors of each element.
        dfdx (:obj:`numpy.array`): Value of the constraint derivative.
        df0dx (:obj:`numpy.array`): Derivative value of the objective function.
        df0dx2 (:obj:`numpy.array`, optional): Derivative value of the second objective function.

    Returns:
        Density filter applied to the derivative values.
    """
    new_deriv_f1, deriv_f, cols = set_deriv(dfdx, df0dx, df0dx2)
    for el in range(deriv_f.shape[0]):
        idx = neighbors[el, :]
        hj  = H1[idx, :].sum(axis=1)
        for ind in range(cols):
            new_deriv_f1[el, ind] = np.sum((1/hj) * H1[el, :] * deriv_f[idx, ind])
    return out_deriv(new_deriv_f1, dfdx, df0dx2)

def normalize(n, f0_scale, f0val, df0dx):
    """ Apply the sensitivity filter to the derivative of the functions (constrain, objective and multiobjective).

    Args:
        n (:obj:`float`): Weight associated with function.
        f0_scale (:obj:`float`): Value of the function.
        f0val (:obj:`numpy.array`): Objective function.
        df0dx (:obj:`numpy.array`): Derivative value of the function.
    
    Returns:
        Normalized function and it's derivative.
    """
    f0val = n * 100 * f0val/f0_scale
    df0dx = n * 100 * df0dx/f0_scale
    return f0val, df0dx

def sensitivity_filter(H, neighbors, x_min_k, xval, dfdx, df0dx, df0dx2=None):
    """ Apply the sensitivity filter to the derivative of the functions (constrain, objective and multiobjective).

    Args:
        H (:obj:`csc_matrix`): Radius subtracted from the distance between the element and the neighbors.
        neighbors (:obj:`csc_matrix`): Neighbors of each element.
        xval (:obj:`numpy.array`): Indicates where there is mass.
        dfdx (:obj:`numpy.array`): Value of the constraint derivative.
        df0dx (:obj:`numpy.array`): Derivative value of the objective function.
        df0dx2 (:obj:`numpy.array`, optional): Derivative value of the second objective function.

    Returns:
        Sensitivity filter applied to the derivative values.
    """
    xval2 = xval.copy()
    xval2[xval2 <= x_min_k] = x_min_k
    
    aux1 = H * xval2[neighbors.flatten()].reshape(H.shape)
    aux3 = 1/(np.sum(H, axis=1) * xval2[:, 0])

    new_deriv_f, deriv_f, cols = set_deriv(dfdx, df0dx, df0dx2)
    for col in range(cols):
        aux2 = aux1 * deriv_f[neighbors.flatten(), col].reshape(H.shape)
        new_deriv_f[:, col] = aux3 * np.sum(aux2, axis=1)
    return out_deriv(new_deriv_f, dfdx, df0dx2)

def apply_constr(fval, dfdx, constr_func, constr_values, freq_constr, ind_constr2, lx, ly, ngl, coord, connect, E, v, rho, alpha_par, beta_par, eta_par, p_par, q_par, x_min_m, x_min_k, area, xval, modes, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, load_vector, omega_par, const_func, free_ind, passive_el, ind_dofs, gradients=True):
    ''' Calculates the constraint functions and derivatives. '''
    if passive_el is not None:
        ind_passive = ind_dofs[passive_el, :]

    for ind in range(len(constr_func)):
        if constr_func[ind] == 'Area':
            aux_fval = total_area(lx, ly, area, xval)
            if gradients:
                aux_dfdx = 100/(lx * ly) * area.reshape(1, -1)

        elif constr_func[ind] == 'Local Ep':
            if ind_constr2[ind] == 0:
                omega_localep = 2 * np.pi * freq_constr[ind]
                dyna_stif_localep = assembly_dyna_stif(omega_localep, mass_matrix, damp_matrix, stif_matrix)
                disp_vector_localep, _, _ = get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif_localep, load_vector, free_ind, omega_localep, alpha_par, beta_par, eta_par, ngl)
            
            aux_fval, fvirg = obj.objective_funcs('Local Ep', disp_vector_localep, passive_el=passive_el, ind_passive=ind_passive, coord=coord, connect=connect, E=E, v=v, rho=rho)
                                                         
            if gradients:
                aux_dfdx = df.derivatives_objective('Local Ep', fvirg, disp_vector_localep, coord, connect, E, v, rho, alpha_par, beta_par, omega_localep, p_par, q_par, x_min_m, x_min_k, xval, \
                           dyna_stif=dyna_stif_localep, ind_dofs=ind_dofs, ngl=ngl, ind_passive=ind_passive, passive_el=passive_el)                     

        elif constr_func[ind] == 'Local Ki':
            if ind_constr2[ind] == 0:
                omega_localki = 2 * np.pi *  freq_constr[ind]
                dyna_stif_localki = assembly_dyna_stif(omega_localki, mass_matrix, damp_matrix, stif_matrix)
                disp_vector_localki, _, _ = get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif_localki, load_vector, free_ind, omega_localki, alpha_par, beta_par, eta_par, ngl)
            
            aux_fval, fvirg = obj.objective_funcs('Local Ki', disp_vector_localki, omega_par=omega_localki, passive_el=passive_el, ind_passive=ind_passive, coord=coord, connect=connect, E=E, v=v, rho=rho)

            if gradients:
                aux_dfdx = df.derivatives_objective('Local Ki', fvirg, disp_vector_localki, coord, connect, E, v, rho, alpha_par, beta_par, omega_localki, p_par, q_par, x_min_m, x_min_k, xval, \
                                dyna_stif=dyna_stif_localki, ind_dofs=ind_dofs, ngl=ngl, ind_passive=ind_passive, passive_el=passive_el)

        elif constr_func[ind] == 'Local R':
            if ind_constr2[ind] == 0:
                omega_localR = 2 * np.pi *  freq_constr[ind]
                dyna_stif_localR = assembly_dyna_stif(omega_localR, mass_matrix, damp_matrix, stif_matrix)
                disp_vector_localR, _, _ = get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif_localR, load_vector, free_ind, omega_localR, alpha_par, beta_par, eta_par, ngl)
            
            aux_fval, fvirg = obj.objective_funcs('Local R', disp_vector_localR, omega_par=omega_localR, passive_el=passive_el, ind_passive=ind_passive, coord=coord, connect=connect, E=E, v=v, rho=rho)

            if gradients:
                aux_dfdx = df.derivatives_objective('Local R', fvirg, disp_vector_localR, coord, connect, E, v, rho, alpha_par, beta_par, omega_localR, p_par, q_par, x_min_m, x_min_k, xval, \
                                dyna_stif=dyna_stif_localR, ind_dofs=ind_dofs, ngl=ngl, ind_passive=ind_passive, passive_el=passive_el)
        
        elif constr_func[ind] == 'Compliance':
            if ind_constr2[ind] == 0:
                omega_comp = 2 * np.pi *  freq_constr[ind]
                dyna_stif_comp = assembly_dyna_stif(omega_comp, mass_matrix, damp_matrix, stif_matrix)
                disp_vector_comp, _, _ = get_disp_vector(modes, stif_matrix, mass_matrix, dyna_stif_comp, load_vector, free_ind, omega_comp, alpha_par, beta_par, eta_par, ngl)

            aux_fval, fvirg = obj.objective_funcs('Compliance', disp_vector_comp, load_vector=load_vector)
            if gradients:
                aux_dfdx = df.derivatives_objective('Compliance', fvirg, disp_vector_comp, coord, connect, E, v, rho, alpha_par, beta_par, omega_comp, p_par, q_par, x_min_m, x_min_k, xval, \
                                                load_vector=load_vector)

        elif constr_func[ind] == 'R Ratio':
            aux_fval, fvirg = obj.objective_funcs('R Ratio', disp_vector, stif_matrix, mass_matrix, load_vector, omega_par, const_func)
            if gradients:
                aux_dfdx = df.derivatives_objective('R Ratio', disp_vector,fvirg, disp_vector, coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, \
                                                mass_matrix=mass_matrix, stif_matrix=stif_matrix, dyna_stif=dyna_stif, free_ind=free_ind)                                                     
        
        if constr_values[ind] < 0:
            aux_fval *= -1
            if gradients:
                aux_dfdx *= -1 
        
        fval[ind, 0] = aux_fval
        if gradients:
            dfdx[ind, :] = aux_dfdx[:, 0]
    return fval, dfdx

def constr_with_freq(constr_func, constr_values):
    ''' Separates the constraint value and the frequency of the constrain function.

    Args:
        constr_values (:obj:`list`): Values of constraint functions applied.
        constr_func (:obj:`list`)  : Constraint functions applied.

    Returns:
        constr_values, freq_constr, ind_freq_constr
    '''
    freq_constr = np.empty(len(constr_func))
    freq_constr.fill(np.nan)
    ind_freq_constr = []
    ind_constr2 = np.zeros(len(constr_func))

    aux_c  = None
    aux_ep = None
    aux_ki = None
    aux_r = None

    for i, value in enumerate(constr_values):
        if type(value) is tuple:
            freq_constr[i] = value[1]
            constr_values[i] = value[0]
            ind_freq_constr.append(i)

            if constr_func[i] == 'Compliance':
                if aux_c is not None and aux_c == value[1]:
                    ind_constr2[i] = 1
                aux_c = value[1]
            elif constr_func[i] == 'Local Ep':
                if aux_ep is not None and aux_ep == value[1]:
                    ind_constr2[i] = 1
                aux_ep = value[1]
            elif constr_func[i] == 'Local Ki':
                if aux_ki is not None and aux_ki == value[1]:
                    ind_constr2[i] = 1
                aux_ki = value[1]
            elif constr_func[i] == 'Local R':
                if aux_r is not None and aux_r == value[1]:
                    ind_constr2[i] = 1
                aux_r = value[1]
    return constr_values, freq_constr, ind_freq_constr, ind_constr2

def normalize_constr(fval, f_scale_constr):
    fval[:, 0] = 100 * fval[:, 0]/f_scale_constr
    return fval

def update_lists(outit, fval, f0val, list_iter, list_fvals, list_f0val, constr_values):
    """ Add new values to list of functions to plot convergence

    Args:
        outit (:obj:`int`): Iteration.
        fval (:obj:`numpy.array`): Constraint function.
        f0val (:obj:`numpy.array`): Objective function.
        list_iter (:obj:`list`): All iteration values.
        list_fvals (:obj:`list`): All constraint function values.
        list_f0val (:obj:`list`): All objective function values.
        constr_values (:obj:`list`): Values of constraint functions applied.

    Returns:
        A tuple of lists with iterations, objective and constraint function.
    """
    list_iter[outit] = outit
    list_f0val[outit] = f0val

    for i in range(len(constr_values)):
        if constr_values[i] > 0:
            list_fvals[outit, i] = (fval[i, 0] + constr_values[i])
        else:
            list_fvals[outit, i] = (fval[i, 0] - constr_values[i])
    return list_iter, list_f0val, list_fvals

def calc_A(coord, ind):
    """ Calculates the total area.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        ind (:obj:`numpy.array`): Element connectivity.

    Returns:
        The element area.
    """   
    area = np.empty((len(ind), 1))
    area[:, 0] = abs(coord[ind[:, 0], 1] * coord[ind[:, 1], 2] + coord[ind[:, 3], 1] * coord[ind[:, 0], 2] + \
            coord[ind[:, 1], 1] * coord[ind[:, 3], 2] - (coord[ind[:, 3], 1] * coord[ind[:, 1], 2] + \
            coord[ind[:, 0], 1] * coord[ind[:, 3], 2] + coord[ind[:, 1], 1] * coord[ind[:, 0], 2]))

    area[:, 0] += abs(coord[ind[:, 0], 1] * coord[ind[:, 3], 2] + coord[ind[:, 2], 1] * coord[ind[:, 0], 2] + \
            coord[ind[:, 3], 1] * coord[ind[:, 2], 2] - (coord[ind[:, 2], 1] * coord[ind[:, 3], 2] + \
            coord[ind[:, 0], 1] * coord[ind[:, 2], 2] + coord[ind[:, 3], 1] * coord[ind[:, 0], 2]))
    area *= 0.5 
    return area
    
def total_area(lx, ly, area, xval):
    """ Calculates the total element area.

    Args:
        lx (:obj:`float`): X-axis length.
        ly (:obj:`float`): Y-axis length.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.

    Returns:
        Total element area.
    """
    return (100/(lx * ly)) * np.sum(xval * area)

def calc_xnew(H, neighbors, xval):
    """ Recalculate xval.

    Args:
        H (:obj:`csc_matrix`): Radius subtracted from the distance between the element and the neighbors.
        neighbors (:obj:`csc_matrix`): Neighbors of each element.
        xval (:obj:`numpy.array`): Indicates where there is mass.
        
    Returns:
        New xval values.
    """
    a = 1/np.sum(H, axis=1)
    b = np.sum(H * xval[neighbors.flatten()].reshape(H.shape), axis=1)
    xe = a * b
    return xe.reshape(-1, 1)

def set_initxval(constr_func, constr_values):
    """ Calculate the initial value of xval.

    Args:
        constr_func (:obj:`list`): constraint functions applied.
        constr_values (:obj:`list`): Values of constraint functions applied.
    
    Returns:
        initial_xval (:obj:`float`): First value of xval
    """
    idx = [i for i, e in enumerate(constr_func) if e == constr_func[0]]
    if len(idx) > 1:
        initial_xval = ((abs(constr_values[idx[0]]) + abs(constr_values[idx[1]]))/2)/100
    else:
        initial_xval = abs(constr_values[idx[0]])/100
    return initial_xval

def get_passive_el(passive_coord, centroids):
    ''' Get index of passive elements .
    
    Args:
        passive_coord (:obj:`tuple`): Region that the shape will not be changed.
        centroids (:obj:`numpy.array`): Coordinate (x,y) of the centroid of each element.

    Returns:
        Index of passive elements.
    '''
    mask = (centroids[:, 0] >= passive_coord[0][0]) & (centroids[:, 0] <= passive_coord[0][1]) & (centroids[:, 1] >= passive_coord[1][0]) & (centroids[:, 1] <= passive_coord[1][1])
    return (mask > 0).nonzero()[0]

def set_passive_el(xmin, xval, passive_el):
    ''' Set the values of passive elements.

    Args:
        xmin (:obj:`numpy.array`): 
        xval (:obj:`numpy.array`): Indicates where there is mass.
        passive_el (:obj:`numpy.array`): Passive element nodes.

    Returns:
        A tuple with updated xmin and xval.
    '''
    if xmin is not None:
        xmin[passive_el] = 0.99
    
    if xval is not None:
        xval[passive_el] = 1
    
    return xmin, xval

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """ Call in a loop to create terminal progress bar.
    https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    
    Args:
        iteration   (:obj:`int`): current iteration 
        total       (:obj:`int`): total iterations 
        prefix      (:obj:`str`, optional): prefix string 
        suffix      (:obj:`str`, optional): suffix string 
        decimals    (:obj:`int`, optional): positive number of decimals in percent complete 
        length      (:obj:`int`, optional): character length of bar 
        fill        (:obj:`str`, optional): bar fill character 
        printEnd    (:obj:`str`, optional): end character (e.g. "\r", "\r\n") 
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

# def finite_difference(mesh_file, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix=None, freq1=180, const_func=100, fac_ratio=2.1, modes=None, rho=7860, E=210e9, v=0.3, x_min_m=0.001, x_min_k=0.001, alpha_par=0, beta_par=5e-6, eta_par=0, p_par=3, q_par=1, passive_coord=None, nodes=[0], number_deltas=5, delta_interval=(1e-12, 1e-2)):
#     ''' Approximate method for solving partial differential equations.
#         Args:
#             el: List of elements to be used.
#                 Example: el = [0, 1, 2] 
#             number_deltas: Number of delta to calculate FDM.
#             delta_interval: Delta range used. 
#     '''
#     l = number_deltas * len(nodes) + 5
#     progress = 0
#     printProgressBar(progress, l, prefix = 'Progress:', suffix = 'Complete', length = 50)

#     low = delta_interval[0]
#     high = delta_interval[1]
#     delta_d = np.linspace(low, high, number_deltas)
       
#     dw      = np.empty((number_deltas, len(nodes)))
#     dw_orig = np.empty(len(nodes))

#     if mesh_file is not None:
#         path = os.path.dirname(os.path.realpath(__file__)) 
#         m_file = os.path.join(path, mesh_file)
#         coord, connect = import_mesh(m_file)
#         ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)
#         nelx = len(coord[coord[:, 2] == coord[0, 2]]) - 1
#         nely = len(coord[coord[:, 1] == coord[0, 1]]) - 1
#         lx = max(coord[:, 1])
#         ly = max(coord[:, 2])
#     else:
#         coord, connect = fc.regularmeshQ4(lx, ly, nelx, nely)
#     ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)
    
#     restri_matrix = fc.get_matrices(restri_matrix, coord, False)
#     free_ind = None
#     if restri_matrix is not None:
#         restricted_ind = fc.get_dofs(restri_matrix)
#         free_ind = fc.remove_dofs(nelx, nely, restricted_ind)
    
#     load_matrix = fc.get_matrices(load_matrix, coord, True)
#     load_vector = fc.get_load_vector(nelx, nely, load_matrix)

#     ngl = 2 * ((nelx + 1) * (nely + 1))
#     xval = 0.5 * np.ones((nelx * nely, 1))
#     if passive_coord is not None:
#         radius = fac_ratio * lx/nelx
#         _, _, centroids = get_neighbors_radius(nelx, nely, coord, connect, radius)
#         passive_el = get_passive_el(passive_coord, centroids)
#         _, xval = set_passive_el(None, xval, passive_el)

#         ind_dofs = fc.get_ind_dofs(connect, 2)
#         ind_passive = ind_dofs[passive_el, :]
#     else:
#         passive_el = None
#         ind_dofs = None
#         ind_passive = None
   
#     omega1_par = 2 * np.pi * freq1  
    
#     data_k, data_m, _ = solution2D(coord, connect, nelx, nely, E, v, rho, xval, x_min_m, x_min_k, p_par, q_par)
#     stif_matrix, mass_matrix, damp_matrix = assembly_matrices(data_k, data_m, ind_rows, ind_cols, ngl, alpha_par, beta_par)
#     dyna_stif = assembly_dyna_stif(omega1_par, mass_matrix, damp_matrix, stif_matrix)
#     if modes is not None:
#         natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)
#         disp_vector, _ = mode_superposition(natural_frequencies, modal_shape, stif_matrix, load_vector, omega1_par, alpha_par, beta_par, eta_par, free_ind)
#     else: 
#         disp_vector, _ = harmonic_problem(ngl, dyna_stif, load_vector, free_ind)
    
#     # F0VAL ORIGINAL
#     f0val_orig, fvirg_orig = objective_funcs(func_name, disp_vector, stif_matrix, mass_matrix, load_vector, omega1_par, const_func, passive_el, ind_passive, coord, connect, E, v, rho)
    
#     # RELACIONADO A DERIVADA    
#     dw_orig[:] = derivatives_objective(func_name, fvirg_orig, disp_vector, coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min_m, x_min_k, xval, \
#                                        load_vector, mass_matrix, stif_matrix, dyna_stif, free_ind, ind_dofs, ngl, passive_el, ind_passive)[nodes, 0]
#     progress += 5
#     printProgressBar(progress, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
#     for i in range(number_deltas):
#         for ind2, node in enumerate(nodes):
#             xval_aux = xval.copy()
#             xval_aux[node] += delta_d[i] 

#             data_k, data_m, _ = solution2D(coord, connect, nelx, nely, E, v, rho, xval_aux, x_min_m, x_min_k, p_par, q_par)
#             stif_matrix, mass_matrix, damp_matrix = assembly_matrices(data_k, data_m, ind_rows, ind_cols, ngl, alpha_par, beta_par)
#             dyna_stif = assembly_dyna_stif(omega1_par, mass_matrix, damp_matrix, stif_matrix)
#             if modes is not None:
#                 natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)
#                 disp_vector, _ = mode_superposition(natural_frequencies, modal_shape, stif_matrix, load_vector, omega1_par, alpha_par, beta_par, eta_par, free_ind)
#             else: 
#                 disp_vector, _ = harmonic_problem(ngl, dyna_stif, load_vector, free_ind)

#             f0val, fvirg = objective_funcs(func_name, disp_vector, stif_matrix, mass_matrix, load_vector, omega1_par, const_func, passive_el, ind_passive, coord, connect, E, v, rho)
            
#             #dw[i, ind2] = (fvirg - fvirg_orig)/delta_d[i]
#             dw[i, ind2] = (f0val - f0val_orig)/delta_d[i]
#             progress += 1
#             printProgressBar(progress, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
#     plt_opt.compare_deriv(nodes, delta_d, dw, dw_orig)

# def freqrsp_multi_modes(modes, freq_rsp, constr_func, constr_values, load_matrix, restri_matrix, coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, x_min_m, x_min_k, p_par, q_par, func_name, const_func, save=False):
#     initial_xval = set_initxval(constr_func, constr_values)
#     xval  = initial_xval * np.ones((nelx * nely, 1))
#     load_vector     = fc.get_load_vector(nelx, nely, load_matrix)
#     restricted_dofs = fc.get_dofs(restri_matrix)
#     free_ind = fc.remove_dofs(nelx, nely, restricted_dofs)
#     ngl = 2 * ((nelx + 1) * (nely + 1))
    
#     freq_range = freq_rsp[:2]
#     delta = freq_rsp[2]    
    
#     rows = int((freq_rsp[1] - freq_rsp[0]) /freq_rsp[2] + 1)
#     func_vector = np.empty((rows, 5), dtype=complex)
#     origin = freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_m, x_min_k, p_par, q_par, freq_range, delta, func_name, const_func, None, load_vector, unrestricted_ind=free_ind)
#     for i, mode in enumerate(modes):      
#         func_vector[:, i] = freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_m, x_min_k, p_par, q_par, freq_range, delta, func_name, const_func, mode, load_vector, unrestricted_ind=free_ind)
#     plt_opt.freqrsp_modes(freq_range, delta, func_vector, origin, modes, func_name, save)

# def freq_test(mesh_file, nelx, nely, lx, ly, E, v, rho, alpha, beta, eta, initial_xval, x_min_m, x_min_k, p_par, q_par, freq_range, func_name, const_func, modes, load_matrix, **kwargs):
#     if mesh_file is not None:
#         path = os.path.dirname(os.path.realpath(__file__)) 
#         m_file = os.path.join(path, mesh_file)
#         coord, connect = import_mesh(m_file)
#         ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)
#         nelx = len(coord[coord[:, 2] == coord[0, 2]]) - 1
#         nely = len(coord[coord[:, 1] == coord[0, 1]]) - 1
#         lx = max(coord[:, 1])
#         ly = max(coord[:, 2])
#     else:
#         coord, connect = fc.regularmeshQ4(lx, ly, nelx, nely)
#     ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)
#     # Force and restrictions matrix
#     load_matrix = fc.get_matrices(load_matrix, coord, True)
#     load_vector = fc.get_load_vector(nelx, nely, load_matrix)
#     if kwargs.get('restri_matrix') is not None: 
#         restri_matrix = kwargs.get('restri_matrix')
#         restri_matrix = fc.get_matrices(restri_matrix, coord, False)
#         restricted_ind = fc.get_dofs(restri_matrix)
#         free_ind = fc.remove_dofs(nelx, nely, restricted_ind)
#     else:
#         free_ind = None
#     ngl = 2 * ((nelx + 1) * (nely + 1))
#     xval = initial_xval * np.ones((nelx * nely, 1))
#     vector = freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_m, x_min_k, p_par, q_par, freq_range[:2], freq_range[2], func_name, const_func, modes, load_vector, unrestricted_ind=free_ind)
#     plt_opt.freqrsp_plot(freq_range, vector)    