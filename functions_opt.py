import cmath
import numpy as np
from time import time
import functions_2d as fc
import plots_opt as plt_opt
from scipy.linalg import eigh
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

def exe_opt(mma, nelx, nely, lx, ly, func_name, load_matrix, restri_matrix=None, freq1=180, constr_func=['Area'], constr_values=[50], n1=1, multiobjective=(None, 0), const_func=100, fac_ratio=2.1, modes=None, rho=7860, E=210e9, v=0.3, x_min_m=0.001, x_min_k=0.001, alpha_par=0, beta_par=5e-6, eta_par=0, alpha_plot=0, beta_plot=1e-8, eta_plot=0, p_par=3, q_par=1, passive_coord=None, freq_rsp=[], chtol=1e-4, dens_filter=True, each_iter=True, max_iter=100, mesh_deform=False, factor=1000, save=False, timing=False):
    if mma:
        import beam_mma as beam
        beam.main(nelx, nely, lx, ly, func_name, load_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, passive_coord, freq_rsp, chtol, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)
    else:
        import beam_gcmma as beam
        beam.main(nelx, nely, lx, ly, func_name, load_matrix, restri_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, passive_coord, freq_rsp, chtol, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)

def solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, omega_par, xval, x_min_m, x_min_k, p_par, q_par):
    """ Assembly matrices.

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
        omega_par (:obj:`float`): 2 pi frequency
        xval (:obj:`numpy.array`): Indicates where there is mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        p_par (:obj:`int`): Penalization power to stiffness.
        q_par (:obj:`int`): Penalization power to mass. 

    Returns:
        A tuple with stiffnes, mass and dynamic matrices and time to assembly the matrices.
    """
    #
    t01 = time() 
    data_k = np.zeros((nelx * nely, 64), dtype=complex)
    data_m = np.zeros((nelx * nely, 64), dtype=complex)
    #
    for el in range(nelx * nely):
        Ke, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        data_k[el, :] = (x_min_k + (xval[el]**p_par)*(1-x_min_k))* Ke.flatten()
        if xval[el]>0.1:
            data_m[el, :] = (x_min_m + (xval[el]**q_par)*(1-x_min_m)) * Me.flatten()
        else:
            data_m[el, :] =  (x_min_m + (3.512e7*xval[el]**9 - 2.081e8*xval[el]**10)*(1-x_min_m) ) * Me.flatten() 
    #
    data_k = data_k.flatten()
    data_m = data_m.flatten()
    stif_matrix = csc_matrix((data_k, (ind_rows, ind_cols)), shape=(ngl, ngl))
    mass_matrix = csc_matrix((data_m, (ind_rows, ind_cols)), shape=(ngl, ngl))
    damp_matrix = alpha * mass_matrix + (beta) * stif_matrix
    dyna_stif = -(omega_par**2) * mass_matrix + 1j * omega_par * damp_matrix + stif_matrix

    tf1 = time()
    t_assembly = str(round((tf1 - t01), 6))
    return stif_matrix.real, mass_matrix.real, dyna_stif, t_assembly

def harmonic_problem(ngl, dyna_stif, load_vector, free_ind=None):
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
        which (:obj:`str`, optional): Which k eigenvectors and eigenvalues to find. 
        sigma (:obj:`float`, optional): Find eigenvalues near sigma using shift-invert mode.
    
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

def mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega_par, alpha, beta, eta, free_ind):    
    """ Perform an harmonic analysis through superposition method and returns the response of
        all nodes due the external or internal equivalent load. It has been implemented two
        different damping models: Viscous Proportional and Hysteretic Proportional
        Entries for Viscous Proportional Model Damping: (alpha_v, beta_v)
        Entries for Hyteretic Proportional Model Damping: (alpha_h, beta_h)
    
    Args:
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        load_vector (:obj:`numpy.array`): Force.
        modes (:obj:`int`, optional): The number of eigenvalues and eigenvectors desired.
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

    natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)

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
    return disp_vector, natural_frequencies, t_superp

def freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_m, x_min_k, p_par, q_par, freq_range, delta, func_name, const_func, modes, load_vector, **kwargs):
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
    free_ind = None
    if kwargs.get('unrestricted_ind') is not None:
        free_ind = kwargs.get('unrestricted_ind')

    interval = np.arange(freq_range[0], freq_range[1] + 1, delta)
    func_vector = np.empty((len(interval)), dtype=complex)

    for n in range(len(interval)):
        omega_par = 2 * np.pi * interval[n]
        stif_matrix, mass_matrix, dyna_stif, _ = solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, omega_par, xval, x_min_m, x_min_k, p_par, q_par)
        if modes is not None:
            disp_vector, _, _ = mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega_par, alpha, beta, eta, free_ind)
        else: 
            disp_vector, _ = harmonic_problem(ngl, dyna_stif, load_vector, free_ind)

        _, func_vector[n] = objective_funcs(func_name, disp_vector, stif_matrix, mass_matrix, load_vector, omega_par, const_func)
        print('It.', n)
    return func_vector 

def compliance(disp_vector, load_vector):
    """ Calculates the compliance function.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        load_vector (:obj:`numpy.array`): Force.
        
    Returns:
        Function value.
    """
    f = abs(np.dot(disp_vector, load_vector))
    return f

def input_power(disp_vector, load_vector, omega_par, const_func):
    """ Calculates the input power function.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        load_vector (:obj:`numpy.array`): Force.
        omega_par (:obj:`float`): 2 * pi * frequency
        const_func (:obj:`float`):

    Returns:
        A tuple with the values of the input power on the logarithmic scale and the 'virgin' input power.
    """
    a = 1j * load_vector.conjugate()@disp_vector
    if omega_par == 0:
        omega_par = 1 #1e-12
    f = 0.5 * omega_par * a.real
    fvirg = f
    #Log Scale
    f = const_func + 10 * np.log10(f.real)
    return f, fvirg

def elastic_potential_energy(disp_vector, stif_matrix, const_func):
    """ Calculates the elastic potential energy.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        const_func (:obj:`float`): 

    Returns:
        A tuple with the values of the potential elastic energy on the logarithmic scale and the 'virgin' potential elastic energy.
    """
    elastic_p = ((1/4) * (disp_vector.reshape(1, -1).conjugate()@stif_matrix@disp_vector))[0]
    fvirg = elastic_p
    #Log Scale
    elastic_p = const_func + 10 * np.log10(elastic_p.real)

    return elastic_p.real, fvirg.real

def kinetic_energy(disp_vector, mass_matrix, omega_par, const_func):
    """ Calculates the kinetic energy.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        omega_par (:obj:`float`): 2 * pi * frequency
        const_func (:obj:`float`):

    Returns:
        A tuple with the values of the kinetic energy on the logarithmic scale and the 'virgin' kinetic energy.
    """
    if omega_par == 0:
        omega_par = 1e-12
    kinetic_e = ((1/4) * omega_par**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
    fvirg = kinetic_e 
    #Log Scale
    kinetic_e  = const_func + 10 * np.log10(kinetic_e)

    return kinetic_e.real, fvirg.real

def R_ratio(disp_vector, stif_matrix, mass_matrix, omega_par, const_func):
    """ Calculates the strain-to-kinetic energy ratio R.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        omega_par (:obj:`float`): 2 * pi * frequency
        const_func (:obj:`float`):

    Returns:
        R, Rvrig 
    """
    elastic_p = ((1/4) * (disp_vector.reshape(1, -1).conjugate()@stif_matrix@disp_vector))[0]
    if omega_par == 0:
        omega_par = 1e-12
    kinetic_e = ((1/4) * omega_par**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
    R = (elastic_p/kinetic_e)
    fvirg = R
    #Log Scale
    R  = const_func + 10 * np.log10(R)
    return R.real, fvirg.real

def lambda_parameter(disp_vector, load_vector, function):
    """ Calculates the lambda parameter of the function.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        load_vector (:obj:`numpy.array`): Force vector.
        function (:obj:`float`): Function value.

    Returns:
        Lambda parameter.
    """
    lam = (disp_vector.conjugate()@load_vector)/function
    return lam

def lambda_parameter_ep(disp_vector, stif_matrix, dyna_stif, free_ind):
    """ Calculates the lambda solution of the elastic potential energy.

    Args:
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix. 
        disp_vector (:obj:`numpy.array`): Displacement.
        free_ind (:obj:`numpy.array`): Free dofs.

    Returns:
        Lambda parameter solution.
    """
    lam = np.zeros(stif_matrix.shape[0], dtype=complex)
    aux = -(1/2) * (stif_matrix[free_ind, :][:, free_ind]@disp_vector[free_ind].conjugate())
    lam[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], aux)
    return lam

def lambda_parameter_ek(disp_vector, mass_matrix, dyna_stif, omega_par, free_ind):
    """ Calculates the lambda solution of the kinetic energy.

    Args:
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        disp_vector (:obj:`numpy.array`): Displacement.
        dyna_stif (array): Stifness matrix.
        omega_par (:obj:`float`): 2 * pi * frequency
        free_ind (:obj:`numpy.array`): Free dofs.

    Returns:
        Lambda parameter solution.
    """
    lam = np.zeros(mass_matrix.shape[0], dtype=complex)
    if omega_par == 0:
        omega_par = 1e-12
    aux = - (omega_par**2) * (mass_matrix[free_ind, :][:, free_ind]@disp_vector[free_ind].conjugate())
    lam[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], aux)
    return lam

def lambda_parameter_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega_par, fvirg, kinetic_e, free_ind):
    """ Calculates the lambda solution of R.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        Kd_matrix (:obj:`numpy.array`): 
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        omega_par (:obj:`float`): 2 * pi * frequency.
        fvirg (:obj:`float`): 'virgin' function value.
        kinetic_e: Kinetic energy.
        free_ind (:obj:`numpy.array`): Free dofs.

    Returns:
        Lambda parameter solution.
    """
    lam = np.zeros(mass_matrix.shape[0], dtype=complex)
    if omega_par == 0:
        omega_par = 1e-12
    aux = - (1/(2*kinetic_e)) * ((stif_matrix[free_ind, :][:, free_ind] - (omega_par**2)*fvirg*mass_matrix[free_ind, :][:, free_ind])@disp_vector[free_ind].conjugate())
    lam[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], aux)
    return lam

def derivative_compliance(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam):
    """ calculates the derivative of the compliance function.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.
        alpha (:obj:`float`): Damping coefficient proportional to mass.
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency
        p_par (:obj:`float`): Penalization power to stiffness. 
        q_par (:obj:`float`): Penalization power to mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness. 
        xval (:obj:`numpy.array`): Indicates where there is mass.
        disp_vector (:obj:`numpy.array`): Displacement.
        lam (:obj:`float`): Lambda parameter.

    Returns:
        Derivative values.
    """
    deriv_f = np.zeros((len(connect), 1))
    dofs = 2
    ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,2]-1, dofs*connect[:,2],
                            dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,4]-1, dofs*connect[:,4]], dtype=int)-1).T
    for el in range(len(connect)):
        Ke, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        ind = ind_dofs[el, :]
        dKe = p_par * (xval[el]**(p_par - 1))*(1-x_min_k) * Ke
        dCe = alpha * Me + beta * dKe
        if xval[el]>0.1:
            dMe = q_par * (xval[el]**(q_par - 1))*(1-x_min_m) * Me
        else:
            dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-x_min_m) ) * Me        
        dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe
        deriv_f[el, 0] = (-lam *(disp_vector[ind].reshape(1, 8)@dKed@disp_vector[ind].reshape(8, 1)))[0,0].real

    return deriv_f

def derivative_input_power(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, fvirg):
    """ calculates the derivative of the input power function.
    
    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.
        alpha (:obj:`float`): Damping coefficient proportional to mass.
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency
        p_par (:obj:`float`): Penalization power to stiffness. 
        q_par (:obj:`float`): Penalization power to mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness. 
        xval (:obj:`numpy.array`): Indicates where there is mass.
        disp_vector (:obj:`numpy.array`): Displacement.
        fvirg (:obj:`float`): Input power function value.

    Returns:
        Derivative values.

    """
    deriv_f = np.zeros((len(connect), 1))
    dofs = 2
    ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,2]-1, dofs*connect[:,2],
                            dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,4]-1, dofs*connect[:,4]], dtype=int)-1).T
    for el in range(len(connect)):
        Ke, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        ind = ind_dofs[el, :]
        dKe = p_par * (xval[el]**(p_par - 1))*(1-x_min_k) * Ke
        dCe = alpha * Me + beta * dKe
        if xval[el]>0.1:
            dMe = q_par * (xval[el]**(q_par - 1))*(1-x_min_m) * Me
        else:
            dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-x_min_m) ) * Me   
        dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe
        a = 1j * (disp_vector[ind].reshape(1, 8)@dKed@disp_vector[ind].reshape(8, 1))[0,0]
        deriv_f[el, 0] = -0.5 * omega_par * a.real
        #Log Scale
        deriv_f[el, 0] = 10.0*deriv_f[el, 0]*np.log10(np.exp(1))/fvirg  
    return deriv_f

def derivative_ep(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam, fvirg):
    """ calculates the derivative of the elastic potential energy function.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.
        alpha (:obj:`float`): Damping coefficient proportional to mass.
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency
        p_par (:obj:`float`): Penalization power to stiffness. 
        q_par (:obj:`float`): Penalization power to mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness. 
        xval (:obj:`numpy.array`): Indicates where there is mass.
        disp_vector (:obj:`numpy.array`): Displacement.
        lam (:obj:`float`): Lambda parameter.
        fvirg (:obj:`float`): Elastic potential energy function value.

    Returns:
        Derivative values.
    """
    deriv_ep = np.zeros((len(connect), 1), dtype=complex)
    dofs = 2
    ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,2]-1, dofs*connect[:,2],
                            dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,4]-1, dofs*connect[:,4]], dtype=int)-1).T
    for el in range(len(connect)):
        Ke, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        ind = ind_dofs[el, :]
        #dKe1 = p_par * (xval[el]**(p_par - 1))*(1-x_min_k) * Ke.conjugate()
        dKe = p_par * (xval[el]**(p_par - 1))*(1-x_min_k) * Ke
        dCe = alpha * Me + beta * dKe
        if xval[el]>0.1:
            dMe = q_par * (xval[el]**(q_par - 1))*(1-x_min_m) * Me
        else:
            dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-x_min_m) ) * Me 
        dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe             
        deriv_ep[el, 0] = (1/4) * (disp_vector[ind].conjugate()@dKe@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        #Log Scale
        deriv_ep[el, 0] = 10.0*deriv_ep[el, 0]*np.log10(np.exp(1))/fvirg
    return deriv_ep.real

def derivative_ek(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam, fvirg):
    """ calculates the derivative of the kinetic energy function.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.
        alpha (:obj:`float`): Damping coefficient proportional to mass.
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency
        p_par (:obj:`float`): Penalization power to stiffness. 
        q_par (:obj:`float`): Penalization power to mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness. 
        xval (:obj:`numpy.array`): Indicates where there is mass.
        disp_vector (:obj:`numpy.array`): Displacement.
        lam (:obj:`float`): Lambda parameter.
        fvirg (:obj:`float`): Kinetic energy function value.

    Returns:
        Derivative values.
    """
    deriv_ek = np.zeros((len(connect), 1), dtype=complex)
    dofs = 2
    ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,2]-1, dofs*connect[:,2],
                            dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,4]-1, dofs*connect[:,4]], dtype=int)-1).T
    for el in range(len(connect)):
        Ke, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        ind = ind_dofs[el, :]
        dKe = p_par * (xval[el]**(p_par - 1))*(1-x_min_k) * Ke
        dCe = alpha * Me + beta * dKe
        if xval[el]>0.1:
            dMe = q_par * (xval[el]**(q_par - 1))*(1-x_min_m) * Me
        else:
            dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-x_min_m) ) * Me 
        dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe             
        deriv_ek[el, 0] = ((omega_par**2)/4) * (disp_vector[ind].conjugate()@dMe@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        #Log Scale
        deriv_ek[el, 0] = 10.0*deriv_ek[el, 0]*np.log10(np.exp(1))/fvirg
    return deriv_ek.real

def derivative_R(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam, fvirg, elastic_p, kinetic_e):
    """ calculates the derivative of the kinetic energy function.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.
        alpha (:obj:`float`): Damping coefficient proportional to mass.
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency
        p_par (:obj:`float`): Penalization power to stiffness. 
        q_par (:obj:`float`): Penalization power to mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness. 
        xval (:obj:`numpy.array`): Indicates where there is mass.
        disp_vector (:obj:`numpy.array`): Displacement.
        lam (:obj:`float`): Lambda parameter.
        fvirg (:obj:`float`): R Ratio function value.
        elastic_p (:obj:`float`): Elastic potential energy function value.
        kinetic_e (:obj:`float`): Kinetic energy function value.

    Returns:
        Derivative values.
    """
    deriv_R = np.zeros((len(connect), 1), dtype=complex)
    dofs = 2
    ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,2]-1, dofs*connect[:,2],
                            dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,4]-1, dofs*connect[:,4]], dtype=int)-1).T
    for el in range(len(connect)):
        Ke, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        ind = ind_dofs[el, :]
        dKe = p_par * (xval[el]**(p_par - 1))*(1-x_min_k) * Ke
        dCe = alpha * Me + beta * dKe
        if xval[el]>0.1:
            dMe = q_par * (xval[el]**(q_par - 1))*(1-x_min_m) * Me
        else:
            dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-x_min_m) ) * Me 
        dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe        
        #
        deriv_R[el, 0] = 1/(4*kinetic_e) * (disp_vector[ind].conjugate()@(dKe - (omega_par**2)*fvirg*dMe)@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        #Log Scale
        deriv_R[el, 0] = 10.0*deriv_R[el, 0]*np.log10(np.exp(1))/fvirg
    return deriv_R.real

def set_deriv(dfdx, df0dx, df0dx2):
    if df0dx2 is not None:
        cols = 2 + dfdx.shape[0]
    else:
        cols = 1 + dfdx.shape[0]
    #
    all_deriv_f = np.empty((df0dx.shape[0], cols))
    for i in range(dfdx.shape[0]):
        all_deriv_f[:, i] = dfdx[i, :]
    #
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

def new_density_filter(H, neighbors, dfdx, df0dx, df0dx2=None):
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
    new_deriv_f, deriv_f, cols = set_deriv(dfdx, df0dx, df0dx2)
    for el in range(deriv_f.shape[0]):
        H_el = H[el, :].data
        idx = neighbors[el, :].data
        Hj = 0
        aux = np.zeros(cols)
        for i in range(len(idx)):
            nn = idx[i]
            Hj = np.sum(H[nn, :])   
            for ind in range(cols):
                aux[ind] +=  (1/Hj) * H_el[i] * deriv_f[nn, ind]
                #aux += (1/Hj) * H_el[i] * deriv_f[nn, 0]
        new_deriv_f[el, :] = aux  
    return out_deriv(new_deriv_f, dfdx, df0dx2)

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

def new_sensitivity_filter(H, neighbors, xval, dfdx, df0dx, df0dx2=None):
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
    aux1 = H.multiply(xval[neighbors.toarray().flatten()].reshape(H.shape))
    aux3 = 1/np.multiply(np.sum(H, axis=1), xval)

    new_deriv_f, deriv_f, cols = set_deriv(dfdx, df0dx, df0dx2)
    for col in range(cols):
        
        aux2 = aux1.multiply(deriv_f[neighbors.toarray().flatten(), col].reshape(H.shape))
        new_deriv_f[:, col] = np.asarray(np.multiply(aux3, np.sum(aux2, axis=1)))[:,0]
    return out_deriv(new_deriv_f, dfdx, df0dx2)

def new_apply_constr(fval, dfdx, constr_func, constr_values, freq_comp_constr, lx, ly,  ind_rows, ind_cols, nelx, nely, coord, connect, E, v, rho, alpha_par, beta_par, eta_par, p_par, q_par, x_min_m, x_min_k, area, xval, modes, disp_vector, dyna_stif, stif_matrix, mass_matrix, load_vector, omega_par, const_func, free_ind, gradients=True):
    ''' Calculates the constraint functions and derivatives. 
    '''
    i = 0
    for ind in range(len(constr_func)):
        if constr_func[ind] == 'Area':
            aux_fval = total_area(lx, ly, area, xval)
            if gradients:
                aux_dfdx = 100/(lx * ly) * area.reshape(1, -1)
        elif constr_func[ind] == 'Compliance':
            if (freq_comp_constr[0] == freq_comp_constr[1]) and (i == 0):
                omega_comp = 2 * np.pi * freq_comp_constr[0]
                ngl = 2 * ((nelx + 1) * (nely + 1))
                stif_matrix_comp, mass_matrix_comp, dyna_stif_comp, _ = solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_par, beta_par, eta_par, omega_comp, xval, x_min_m, x_min_k, p_par, q_par)
                if modes is not None:
                    disp_vector_comp, _, _ = mode_superposition(stif_matrix_comp, mass_matrix_comp, load_vector, modes, omega_comp, alpha_par, beta_par, eta_par, free_ind)
                else: 
                    disp_vector_comp, _ = harmonic_problem(ngl, dyna_stif_comp, load_vector, free_ind)
            elif not(freq_comp_constr[0] == freq_comp_constr[1]):    
                omega_comp = 2 * np.pi * freq_comp_constr[i]
                ngl = 2 * ((nelx + 1) * (nely + 1))
                stif_matrix_comp, mass_matrix_comp, dyna_stif_comp, _ = solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha_par, beta_par, eta_par, omega_comp, xval, x_min_m, x_min_k, p_par, q_par)
                if modes is not None:
                    disp_vector_comp, _, _ = mode_superposition(stif_matrix_comp, mass_matrix_comp, load_vector, modes, omega_comp, alpha_par, beta_par, eta_par, free_ind)
                else: 
                    disp_vector_comp, _ = harmonic_problem(ngl, dyna_stif_comp, load_vector, free_ind)
            aux_fval, fvirg = objective_funcs('Compliance', disp_vector_comp, stif_matrix_comp, mass_matrix_comp, load_vector, omega_comp, const_func)
            if gradients:
                aux_dfdx = derivatives_objective('Compliance', disp_vector_comp, stif_matrix_comp, dyna_stif_comp, mass_matrix_comp, load_vector, fvirg, coord, connect, \
                                                    E, v, rho, alpha_par, beta_par, omega_comp, p_par, q_par, x_min_m, x_min_k, xval, free_ind)
            i += 1
        elif constr_func[ind] == 'R Ratio':
            aux_fval, fvirg = objective_funcs('R Ratio', disp_vector, stif_matrix, mass_matrix, load_vector, omega_par, const_func)
            if gradients:
                aux_dfdx = derivatives_objective('R Ratio', disp_vector, stif_matrix, dyna_stif, mass_matrix, load_vector, fvirg, coord, connect, \
                                                    E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, free_ind)
        #
        if constr_values[ind] < 0:
            aux_fval *= -1
            if gradients:
                aux_dfdx *= -1 
        aux_fval -= constr_values[ind]
        #
        fval[ind, 0] = aux_fval
        if gradients:
            dfdx[ind, :] = aux_dfdx[:, 0]
    return fval, dfdx

def constr_compliance(constr_values, constr_func):
    ''' Separates the constraint value and the frequency of the compliance constraint.
    Args:
        constr_values (:obj:`list`): Values of constraint functions applied.
        constr_func (:obj:`list`)  : Constraint functions applied.

    Returns:
        constr_values, freq_comp_constr, ind_comp
    '''
    freq_comp_constr    = np.empty(2)
    freq_comp_constr[:] = None
    ind_comp = []
    ind      = 0
    for i in range(len(constr_func)):
        if constr_func[i] == 'Compliance':
            freq_comp_constr[ind] = constr_values[i][1]
            constr_values[i]      = constr_values[i][0]
            ind_comp.append(i)
            ind += 1
    return constr_values, freq_comp_constr, ind_comp

def update_lists(outit, fval, f0val, list_iter, list_fvals, list_f0val, constr_func, constr_values):
    """ Add new values to list of functions to plot convergence

    Args:
        outit (:obj:`int`): Iteration.
        fval (:obj:`numpy.array`): Constraint function.
        f0val (:obj:`numpy.array`): Objective function.
        list_iter (:obj:`list`): All iteration values.
        list_fvals (:obj:`list`): All constraint function values.
        list_f0val (:obj:`list`): All objective function values.
        constr_func (:obj:`list`): constraint functions applied.
        constr_values (:obj:`list`): Values of constraint functions applied.

    Returns:
        A tuple of lists with iterations, objective and constraint function.
    """
    list_iter.append(outit)
    list_f0val.append(f0val)
    for i in range(len(constr_func)):
        if constr_values[i] > 0:
            list_fvals[i].append(fval[i, 0] + constr_values[i])
        else:
            list_fvals[i].append(fval[i, 0] - constr_values[i])
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
    # 
    centroids = np.empty((el_number, 2))
    idx = connect[:, 1:] - 1   
    centroids[:, 0] = np.sum(coord[idx, 1], axis = 1)/4
    centroids[:, 1] = np.sum(coord[idx, 2], axis = 1)/4
    # 
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
        #
        hi      = radius - distance
        hi_max  = np.maximum(0, hi)
        data.extend(hi_max[mask])
        aux     = len(hi_max[mask])
        rows    = np.repeat(el, aux) #.tolist()
        columns = np.arange(0, aux)
        ind_rows.extend(rows) 
        ind_cols.extend(columns)
        #
        if aux > cols:
            cols = aux
    H = csc_matrix((data, (ind_rows, ind_cols)), shape=(nelx*nely, cols))
    neighbors = csc_matrix((neighbors, (ind_rows, ind_cols)), shape=(nelx*nely, cols), dtype='int')
    return neighbors, H, centroids

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
    b = np.sum(H.multiply(xval[neighbors.toarray().flatten()].reshape(H.shape)), axis=1)
    xe = np.multiply(a, b)
    return np.asarray(xe)

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

def set_passive_el(xmin, xval, passive_coord, centroids):
    ''' Set the values of passive elements.
        Args:
            xmin (:obj:`numpy.array`): 
            xval (:obj:`numpy.array`): Indicates where there is mass.
            passive_coord (:obj:`tuple`): Region that the shape will not be changed.
                Example: ((0.5, 1), (0.3, 0.6)) = ((x_initial, x_final), (y_initial, y_final))
            centroids (:obj:`numpy.array`): Coordinate (x,y) of the centroid of each element.

        Returns:
            A tuple with updated xmin and xval.
    '''
    mask = (centroids[:, 0] >= passive_coord[0][0]) & (centroids[:, 0] <= passive_coord[0][1]) & (centroids[:, 1] >= passive_coord[1][0]) & (centroids[:, 1] <= passive_coord[1][1])
    xmin[mask] = 0.99
    xval[mask] = 1
    return xmin, xval

def objective_funcs(func_name, disp_vector, stif_matrix, mass_matrix, load_vector, omega_par, const_func):
    ''' Calculate objective function.

    Args:
        func_name (:obj:`str`): Objective function used.
        disp_vector (:obj:`numpy.array`): Displacement.
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        load_vector (:obj:`numpy.array`): Force.
        omega_par (:obj:`float`): 2 * pi * frequency.
        const_func (:obj:`float`):

    Returns:
        Objective function values.
    '''
    if func_name == "Compliance":
        f0val = compliance(disp_vector, load_vector)
        fvirg = f0val
    #
    elif func_name == "Elastic Potential Energy":
        f0val, fvirg = elastic_potential_energy(disp_vector, stif_matrix, const_func)
    #
    elif func_name == "Input Power":
        f0val, fvirg = input_power(disp_vector, load_vector, omega_par, const_func)
    #               
    elif func_name == "Kinetic Energy":
        f0val, fvirg = kinetic_energy(disp_vector, mass_matrix, omega_par, const_func)
    #   
    elif func_name == 'R Ratio':
        f0val, fvirg = R_ratio(disp_vector, stif_matrix, mass_matrix, omega_par, const_func)
    return f0val, fvirg
       
def derivatives_objective(func_name, disp_vector, stif_matrix, dyna_stif, mass_matrix, load_vector, fvirg, coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xnew, free_ind):
    ''' Calculate derivatite of objective function.

    Args:
        func_name (:obj:`str`): Objective function used.
        disp_vector (:obj:`numpy.array`): Displacement.
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        load_vector (:obj:`numpy.array`): Force.
        fvirg (:obj:`float`): 'virgin' function value.
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio. 
        rho (:obj:`float`): Density.
        alpha_par (:obj:`float`): Damping coefficient proportional to mass.
        beta_par (:obj:`float`): Damping coefficient proportional to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency.  
        p_par (:obj:`int`): Penalization power to stiffness.
        q_par (:obj:`int`): Penalization power to mass.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness. 
        xnew (:obj:`numpy.array`): Indicates where there is mass.
        free_ind (:obj:`numpy.array`): DOFs free.

    Returns:
        Derivative values.
    '''
    if func_name == "Compliance":
        lam_par = lambda_parameter(disp_vector, load_vector, fvirg)
        df0dx   = derivative_compliance(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xnew, disp_vector, lam_par)
    #
    elif func_name == "Elastic Potential Energy":
        lam_par = lambda_parameter_ep(disp_vector, stif_matrix, dyna_stif, free_ind)
        df0dx   = derivative_ep(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xnew, disp_vector, lam_par, fvirg)
    #
    elif func_name == "Input Power":
        df0dx = derivative_input_power(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xnew, disp_vector, fvirg)
    #               
    elif func_name == "Kinetic Energy":
        lam_par = lambda_parameter_ek(disp_vector, mass_matrix, dyna_stif, omega_par, free_ind)
        df0dx   = derivative_ek(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xnew, disp_vector, lam_par, fvirg)
    #   
    elif func_name == 'R Ratio':
        elastic_p = ((1/4) * (disp_vector.reshape(1, -1).conjugate()@stif_matrix@disp_vector))[0]
        if omega_par == 0:
            omega_par = 1e-12
        kinetic_e = ((1/4) * omega_par**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
        lam_par   = lambda_parameter_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega_par, fvirg, kinetic_e, free_ind)
        df0dx     = derivative_R(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xnew, disp_vector, lam_par, fvirg, elastic_p, kinetic_e)
    return df0dx

def freqrsp_multi_modes(modes, freq_rsp, constr_func, constr_values, load_matrix, restri_matrix, coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, x_min_m, x_min_k, p_par, q_par, func_name, const_func, save=False):

    initial_xval = set_initxval(constr_func, constr_values)
    xval  = initial_xval * np.ones((nelx * nely, 1))
    load_vector     = fc.get_load_vector(nelx, nely, load_matrix)
    restricted_dofs = fc.get_dofs(restri_matrix)
    free_ind = fc.remove_dofs(nelx, nely, restricted_dofs)
    ngl = 2 * ((nelx + 1) * (nely + 1))
    #
    freq_range = freq_rsp[:2]
    delta = freq_rsp[2]    
    #
    rows = int((freq_rsp[1] - freq_rsp[0]) /freq_rsp[2] + 1)
    func_vector = np.empty((rows, 5), dtype=complex)
    origin = freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_m, x_min_k, p_par, q_par, freq_range, delta, func_name, const_func, None, load_vector, unrestricted_ind=free_ind)
    for i, mode in enumerate(modes):      
        func_vector[:, i] = freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_m, x_min_k, p_par, q_par, freq_range, delta, func_name, const_func, mode, load_vector, unrestricted_ind=free_ind)
    plt_opt.freqrsp_modes(freq_range, delta, func_vector, origin, modes, func_name, save)
    