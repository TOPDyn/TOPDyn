from time import time
import numpy as np
import cmath
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigs, spsolve
import sys
import os
sys.path.insert(1, os.getcwd() + '\SolverFEM2D')
import functions2D as fc

def solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, omega, xval, x_min_k, x_min_m, p_par, q_par):
    """ Assembly matrices.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        ind_rows (numpy.array): Node indexes to make the assembly.
        ind_cols (numpy.array): Node indexes to make the assembly.
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        ngl (int): Degrees of freedom.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.  
        rho (float): Density.  
        alpha (float): Damping coefficient proportional to mass. 
        beta (float): Damping coefficient proportional to stiffness.  
        eta (float): Damping coefficient. 
        omega (float): 2 pi frequency
        xval (numpy.array): Indicates where there is mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        p_par (int): Penalization power to stiffness.
        q_par (int): Penalization power to mass. 

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
    dyna_stif = -(omega**2) * mass_matrix + 1j * omega * damp_matrix + stif_matrix

    tf1 = time()
    t_assembly = str(round((tf1 - t01), 6))
    
    return stif_matrix, mass_matrix, dyna_stif, t_assembly

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

def modal_analysis(stif_matrix, mass_matrix, modes=20, which='LM', sigma=0.01):
    """ Modal Analysis. Use eigs Scipy function.

    Args:
        stif_matrix (numpy.array): Stiffness matrix.
        mass_matrix (numpy.array): Mass matrix.
        modes (:obj:`int`, optional): The number of eigenvalues and eigenvectors desired.
        which (:obj:`str`, optional): Which k eigenvectors and eigenvalues to find. 
        sigma (:obj:`float`, optional): Find eigenvalues near sigma using shift-invert mode.
    
    Returns:
        A tuple with natural frequencies and modes shape.
    """
    eigen_values, eigen_vectors = eigs(stif_matrix, M=mass_matrix, k=modes, which=which, sigma=sigma)

    positive_real = np.absolute(np.real(eigen_values))
    natural_frequencies = np.sqrt(positive_real)/(2*np.pi)
    modal_shape = np.real(eigen_vectors)

    index_order = np.argsort(natural_frequencies)
    natural_frequencies = natural_frequencies[index_order]
    modal_shape = modal_shape[:, index_order]

    return natural_frequencies, modal_shape

def mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega, alpha, beta, eta, free_ind):    
    """ Perform an harmonic analysis through superposition method and returns the response of
        all nodes due the external or internal equivalent load. It has been implemented two
        different damping models: Viscous Proportional and Hysteretic Proportional
        Entries for Viscous Proportional Model Damping: (alpha_v, beta_v)
        Entries for Hyteretic Proportional Model Damping: (alpha_h, beta_h)
    
    Args:
        stif_matrix (numpy.array): Stiffness matrix.
        mass_matrix (numpy.array): Mass matrix.
        load_vector (numpy.array): Force.
        modes (:obj:`int`, optional): The number of eigenvalues and eigenvectors desired.
        omega (float): 2 pi frequency
        alpha (float): Damping coefficient proportional to mass. 
        beta (float): Damping coefficient proportional to stiffness.
        eta (float): Damping coefficient. 
        free_ind (numpy.array): Free dofs. 

    Returns:
        A tuple with displacement and time to solve the problem.
    """
    t0 = time()
    alphaV, betaV, betaH = alpha, beta, eta

    natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)

    F_aux = modal_shape.T @ load_vector[free_ind]
    omega_n = 2*np.pi*natural_frequencies
    F_kg = (omega_n**2)

    F_mg =  - (omega**2)
    F_cg = 1j*((betaH + betaV*omega)*(omega_n**2) + (0. + omega*alphaV)) 
    data = np.divide(1, (F_kg + F_mg + F_cg))
    diag = np.diag(data)
    #disp_vector = modal_shape @ (diag @ F_aux[:,i])
    rows = stif_matrix.shape[0]
    disp_vector = np.zeros((rows), dtype=complex)
    disp_vector[free_ind] = modal_shape @ (diag @ F_aux)

    tf = time()
    t_superp = str(round((tf - t0), 6))
    
    return disp_vector, t_superp

def freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, xval, x_min_k, x_min_m, p_par, q_par, freq_range, delta, function_name, const_func, modes, load_vector, **kwargs):
    """ Calculates the objective function for a range of frequencies.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        ind_rows (numpy.array): Node indexes to make the assembly.
        ind_cols (numpy.array): Node indexes to make the assembly.
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        ngl (int): Degrees of freedom.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.  
        rho (float): Density.  
        alpha (float): Damping coefficient proportional to mass. 
        beta (float): Damping coefficient proportional to stiffness.  
        eta (float): Damping coefficient. 
        xval (numpy.array): Indicates where there is mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        p_par (int): Penalization power to stiffness. 
        q_par (int): Penalization power to mass.
        freq_range (list): Frequency range.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (int): Step between each calculation of the objective function. 
        func_name (str): Objective function used.
        const_func ()
        modes (int): The number of eigenvalues and eigenvectors desired.
        load_vector (numpy.array): Force.

    Returns:
        Objective function values.
    """
    free_ind = None
    if kwargs.get('unrestricted_ind') is not None:
        free_ind = kwargs.get('unrestricted_ind')

    interval = np.arange(freq_range[0], freq_range[1] + 1, delta)
    func_vector = np.empty((len(interval)), dtype=complex)

    for n in range(len(interval)):
        omega = 2 * np.pi * interval[n]
        stif_matrix, mass_matrix, dyna_stif, _ = solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, omega, xval, x_min_k, x_min_m, p_par, q_par)
        if modes is not None:
            disp_vector, _ = mode_superposition(stif_matrix, mass_matrix, load_vector, modes, omega, alpha, beta, eta, free_ind)
        else: 
            disp_vector, _ = harmonic_problem(ngl, dyna_stif, load_vector, free_ind)
        if function_name == "Compliance":
            abs(disp_vector@load_vector)
        elif function_name == "Elastic Potential Energy":
            _, func_vector[n] = elastic_potential_energy(disp_vector, stif_matrix, const_func) 
        elif function_name == "Input Power":
            _, func_vector[n] = func_input_power(disp_vector, load_vector, omega, const_func)
        elif function_name == "Kinetic Energy":
            _, func_vector[n] = kinetic_energy(disp_vector, mass_matrix, omega, const_func)
        elif function_name == 'R Ratio':
            _, func_vector[n], _, _ = R_ratio(disp_vector, stif_matrix, mass_matrix, omega, const_func) 

        print('It.', n)
    return func_vector 

def func_compliance(disp_vector, load_vector):
    """ Calculates the compliance function.

    Args:
        disp_vector (numpy.array): Displacement.
        load_vector (numpy.array): Force.
        
    Returns:
        Function value.
    """
    f = abs(np.dot(disp_vector, load_vector))
 
    return f

def func_input_power(disp_vector, load_vector, omega, const_func):
    """ Calculates the input power function.

    Args:
        disp_vector (numpy.array): Displacement.
        load_vector (numpy.array): Force.
        omega (float): 2 * pi * frequency
        const_func (float):

    Returns:
        A tuple with the values of the input power on the logarithmic scale and the 'virgin' input power.
    """
    a = 1j * load_vector.conjugate()@disp_vector
    if omega == 0:
        omega = 1 #1e-12
    f = 0.5 * omega * a.real
    fvirg = f
    #Log Scale
    f = const_func + 10 * np.log10(f.real)
    return f, fvirg

def elastic_potential_energy(disp_vector, stif_matrix, const_func):
    """ Calculates the elastic potential energy.

    Args:
        disp_vector (numpy.array): Displacement.
        stif_matrix (numpy.array): Stiffness matrix.
        const_func (float): 

    Returns:
        A tuple with the values of the potential elastic energy on the logarithmic scale and the 'virgin' potential elastic energy.
    """
    elastic_p = ((1/4) * (disp_vector.reshape(1, -1).conjugate()@stif_matrix@disp_vector))[0]
    fvirg = elastic_p
    #Log Scale
    elastic_p = const_func + 10 * np.log10(elastic_p.real)

    return elastic_p.real, fvirg.real

def kinetic_energy(disp_vector, mass_matrix, omega, const_func):
    """ Calculates the kinetic energy.

    Args:
        disp_vector (numpy.array): Displacement.
        mass_matrix (numpy.array): Mass matrix.
        omega (float): 2 * pi * frequency
        const_func (float):

    Returns:
        A tuple with the values of the kinetic energy on the logarithmic scale and the 'virgin' kinetic energy.
    """
    if omega == 0:
        omega = 1e-12
    kinetic_e = ((1/4) * omega**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
    fvirg = kinetic_e 
    #Log Scale
    kinetic_e  = const_func + 10 * np.log10(kinetic_e)

    return kinetic_e.real, fvirg.real

def R_ratio(disp_vector, stif_matrix, mass_matrix, omega, const_func):
    """ Calculates the strain-to-kinetic energy ratio R.

    Args:
        disp_vector (numpy.array): Displacement.
        stif_matrix (numpy.array): Stiffness matrix.
        mass_matrix (numpy.array): Mass matrix.
        omega (float): 2 * pi * frequency
        const_func (float):

    Returns:
        R, Rvrig Ep and Ek
    """
    elastic_p = ((1/4) * (disp_vector.reshape(1, -1).conjugate()@stif_matrix@disp_vector).real)[0]
    if omega == 0:
        omega = 1e-12
    kinetic_e = ((1/4) * (omega**2) * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
    R = (elastic_p/kinetic_e)
    fvirg = R
    #Log Scale
    R  = const_func + 10 * np.log10(R)

    return R, fvirg.real, elastic_p.real, kinetic_e.real

def lambda_parameter(disp_vector, load_vector, function):
    """ Calculates the lambda parameter of the function.

    Args:
        disp_vector (numpy.array): Displacement.
        load_vector (numpy.array): Force vector.
        function (float): Function value.

    Returns:
        Lambda parameter.
    """
    lam = (disp_vector.conjugate()@load_vector)/function
    return lam

def lambda_parameter_ep(disp_vector, stif_matrix, dyna_stif, free_ind):
    """ Calculates the lambda solution of the elastic potential energy.

    Args:
        stif_matrix (numpy.array): Stiffness matrix.
        dyna_stif (numpy.array): Dynamic stiffness matrix. 
        disp_vector (numpy.array): Displacement.
        free_ind (numpy.array): Free dofs.

    Returns:
        Lambda parameter solution.
    """
    lam = np.zeros(stif_matrix.shape[0], dtype=complex)
    
    aux = -(1/2) * (stif_matrix[free_ind, :][:, free_ind]@disp_vector[free_ind].conjugate())

    lam[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], aux)

    return lam

def lambda_parameter_ek(disp_vector, mass_matrix, dyna_stif, omega, free_ind):
    """ Calculates the lambda solution of the kinetic energy.

    Args:
        mass_matrix (numpy.array): Mass matrix.
        disp_vector (numpy.array): Displacement.
        dyna_stif (array): Stifness matrix.
        omega (float): 2 * pi * frequency
        free_ind (numpy.array): Free dofs.

    Returns:
        Lambda parameter solution.
    """
    lam = np.zeros(mass_matrix.shape[0], dtype=complex)
    if omega == 0:
        omega = 1e-12
    aux = - (omega**2) * (mass_matrix[free_ind, :][:, free_ind]@disp_vector[free_ind].conjugate())

    lam[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], aux)

    return lam

def lambda_parameter_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega, fvirg, kinetic_e, free_ind):
    """ Calculates the lambda solution of R.

    Args:
        disp_vector (numpy.array): Displacement.
        Kd_matrix (numpy.array): 
        stif_matrix (numpy.array): Stiffness matrix.
        mass_matrix (numpy.array): Mass matrix.
        omega (float): 2 * pi * frequency.
        fvirg (float): 'virgin' function value.
        kinetic_e: Kinetic energy.
        free_ind (numpy.array): Free dofs.

    Returns:
        Lambda parameter solution.
    """
    lam = np.zeros(mass_matrix.shape[0], dtype=complex)
    if omega == 0:
        omega = 1e-12

    aux = - (1/(2*kinetic_e)) * ((stif_matrix[free_ind, :][:, free_ind] - (omega**2)*fvirg*mass_matrix[free_ind, :][:, free_ind])@disp_vector[free_ind].conjugate())
   
    lam[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], aux)

    return lam

def derivative_compliance(coord, connect, E, v, rho, alpha, beta, omega, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, lam):
    """ calculates the derivative of the compliance function.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.
        rho (float): Density.
        alpha (float): Damping coefficient proportional to mass.
        beta (float): Damping coefficient proportional to stiffness.
        omega (float): 2 * pi * frequency
        p_par (float): Penalization power to stiffness. 
        q_par (float): Penalization power to mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        xval (numpy.array): Indicates where there is mass.
        disp_vector (numpy.array): Displacement.
        lam (complex): Lambda parameter.

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
        dKed = dKe + omega * 1j * dCe - (omega**2) * dMe
        deriv_f[el, 0] = (-lam *(disp_vector[ind].reshape(1, 8)@dKed@disp_vector[ind].reshape(8, 1)))[0,0].real

    return deriv_f

def derivative_input_power(coord, connect, E, v, rho, alpha, beta, omega, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, fvirg):
    """ calculates the derivative of the input power function.
    
    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.
        rho (float): Density.
        alpha (float): Damping coefficient proportional to mass.
        beta (float): Damping coefficient proportional to stiffness.
        omega (float): 2 * pi * frequency
        p_par (float): Penalization power to stiffness. 
        q_par (float): Penalization power to mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        xval (numpy.array): Indicates where there is mass.
        disp_vector (numpy.array): Displacement.
        fvirg (float): Input power function value.

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
        dKed = dKe + omega * 1j * dCe - (omega**2) * dMe
        a = 1j * (disp_vector[ind].reshape(1, 8)@dKed@disp_vector[ind].reshape(8, 1))[0,0]
        deriv_f[el, 0] = -0.5 * omega * a.real
        #Log Scale
        deriv_f[el, 0] = 10.0*deriv_f[el, 0]*np.log10(np.exp(1))/fvirg  
    return deriv_f

def derivative_ep(coord, connect, E, v, rho, alpha, beta, omega, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, lam, fvirg):
    """ calculates the derivative of the elastic potential energy function.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.
        rho (float): Density.
        alpha (float): Damping coefficient proportional to mass.
        beta (float): Damping coefficient proportional to stiffness.
        omega (float): 2 * pi * frequency
        p_par (float): Penalization power to stiffness. 
        q_par (float): Penalization power to mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        xval (numpy.array): Indicates where there is mass.
        disp_vector (numpy.array): Displacement.
        lam (complex): Lambda parameter.
        fvirg (float): Elastic potential energy function value.

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
        dKed = dKe + omega * 1j * dCe - (omega**2) * dMe             
        deriv_ep[el, 0] = (1/4) * (disp_vector[ind].conjugate()@dKe@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        #Log Scale
        deriv_ep[el, 0] = 10.0*deriv_ep[el, 0]*np.log10(np.exp(1))/fvirg

    return deriv_ep.real

def derivative_ek(coord, connect, E, v, rho, alpha, beta, omega, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, lam, fvirg):
    """ calculates the derivative of the kinetic energy function.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.
        rho (float): Density.
        alpha (float): Damping coefficient proportional to mass.
        beta (float): Damping coefficient proportional to stiffness.
        omega (float): 2 * pi * frequency
        p_par (float): Penalization power to stiffness. 
        q_par (float): Penalization power to mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        xval (numpy.array): Indicates where there is mass.
        disp_vector (numpy.array): Displacement.
        lam (complex): Lambda parameter.
        fvirg (float): Kinetic energy function value.

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
        dKed = dKe + omega * 1j * dCe - (omega**2) * dMe             
        deriv_ek[el, 0] = ((omega**2)/4) * (disp_vector[ind].conjugate()@dMe@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        #Log Scale
        deriv_ek[el, 0] = 10.0*deriv_ek[el, 0]*np.log10(np.exp(1))/fvirg

    return deriv_ek.real

def derivative_R(coord, connect, E, v, rho, alpha, beta, omega, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, lam, fvirg, elastic_p, kinetic_e):
    """ calculates the derivative of the kinetic energy function.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.
        rho (float): Density.
        alpha (float): Damping coefficient proportional to mass.
        beta (float): Damping coefficient proportional to stiffness.
        omega (float): 2 * pi * frequency
        p_par (float): Penalization power to stiffness. 
        q_par (float): Penalization power to mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        xval (numpy.array): Indicates where there is mass.
        disp_vector (numpy.array): Displacement.
        lam (complex): Lambda parameter.
        fvirg (float): R Ratio function value.
        elastic_p (float): Elastic potential energy function value.
        kinetic_e (float): Kinetic energy function value.

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
        dKed = dKe + omega * 1j * dCe - (omega**2) * dMe        
        #
        deriv_R[el, 0] = 1/(4*kinetic_e) * (disp_vector[ind].conjugate()@(dKe - (omega**2)*fvirg*dMe)@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        #Log Scale
        deriv_R[el, 0] = 10.0*deriv_R[el, 0]*np.log10(np.exp(1))/fvirg

    return deriv_R.real

def dfAt_density(deriv_At, H, neighbors):
    """ Apply the density filter to the derivative of the area.

    Args:
        deriv_At (numpy.array): Derivative of total area.
        H (csc_matrix): Radius subtracted from the distance between the element and the neighbors.
        neighbors (csc_matrix): Neighbors of each element.

    Returns:
        Density filter applied to the derivative values.
    """
    new_deriv_At = np.empty(deriv_At.shape)
    for el in range(len(deriv_At)):
        H_el = H[el, :].data
        idx = neighbors[el, :].data
        Hj = 0
        aux = 0
        for i in range(len(idx)):
            nn = idx[i]
            Hj = np.sum(H[nn, :])   
            aux += (1/Hj) * H_el[i] * deriv_At[nn]
        new_deriv_At[el] = aux

    return new_deriv_At

def density_filter(deriv_f, H, neighbors):
    """ Apply the density filter to the derivative of the function.

    Args:
        deriv_f (numpy.array): Derivative of the function.
        H (csc_matrix): Radius subtracted from the distance between the element and the neighbors.
        neighbors (csc_matrix): Neighbors of each element.

    Returns:
        Density filter applied to the derivative values.
    """
    new_deriv_f = np.empty(deriv_f.shape)    
    for el in range(deriv_f.shape[0]):
        H_el = H[el, :].data
        idx = neighbors[el, :].data
        Hj = 0
        aux = 0
        for i in range(len(idx)):
            nn = idx[i]
            Hj = np.sum(H[nn, :])   
            aux += (1/Hj) * H_el[i] * deriv_f[nn, 0]
        new_deriv_f[el] = aux
          
    return new_deriv_f

def sensitivity_filter(deriv_f, H, neighbors, xval, radius):
    """ Apply the sensitivity filter to the derivative of the function.

    Args:
        deriv_f (numpy.array): Derivative of the function.
        H (csc_matrix): Radius subtracted from the distance between the element and the neighbors.
        neighbors (csc_matrix): Neighbors of each element.
        xval (numpy.array): Indicates where there is mass.
        radius (float): Radius to get elements in the vicinity of each element. 

    Returns:
        Sensitivity filter applied to the derivative values.
    """
    aux1 = H.multiply(xval[neighbors.toarray().flatten()].reshape(H.shape))
    aux2 = aux1.multiply(deriv_f[neighbors.toarray().flatten()].reshape(H.shape))

    aux3 = 1/np.multiply(np.sum(H, axis=1), xval)
    new_deriv = np.multiply(aux3, np.sum(aux2, axis=1))

    return np.asarray(new_deriv)

def dens_dconstr(dfdx, constr_func, H, neighbors, radius):
    """ Apply the density filter to the derivative of the constraint.

    Args:
        dfdx (numpy.array): Derivative of the constraint.
        constr_func (list): Restriction functions applied.
        H (csc_matrix): Radius subtracted from the distance between the element and the neighbors.
        neighbors (csc_matrix): Neighbors of each element.
        xval (numpy.array): Indicates where there is mass.
        radius (float): Radius to get elements in the vicinity of each element.

    Returns:
        Density filter applied to the derivative values.
    """
    for i in range(len(constr_func)):
        if constr_func[i] == 'Area':
            dfdx[i, :] = dfAt_density(dfdx[i, :], H, neighbors)

        if constr_func[i] == 'R Ratio':
            dfdx[i, :] = density_filter(dfdx[i, :].reshape(-1, 1), H, neighbors).reshape(-1)            

    return dfdx

def sens_dconst(dfdx, constr_func, H, neighbors, xval, radius):
    """ Apply the sensitivity filter to the derivative of the constraint.

    Args:
        dfdx (numpy.array): Derivative of the constraint.
        constr_func (list): Restriction functions applied.
        H (csc_matrix): Radius subtracted from the distance between the element and the neighbors.
        neighbors (csc_matrix): Neighbors of each element.
        xval (numpy.array): Indicates where there is mass.
        radius (float): Radius to get elements in the vicinity of each element.

    Returns:
        Sensitivity filter applied to the derivative values.
    """
    for i in range(len(constr_func)):
        if constr_func[i] == 'R Ratio':
            dfdx[i, :] = sensitivity_filter(dfdx[i, :], H, neighbors, xval, radius)[:, 0]

    return dfdx

def area_constr(fval, dfdx, ind, constr_values, lx, ly, area, xval):
    """ Calculates the function and derivative of the area.

    Args:
        fval (numpy.array): Value of the constraint function.
        dfdx (numpy.array): Value of the constraint derivative.
        ind (int): Function index in the constr_func list.
        constr_values (list): Values of restriction functions applied. 
        lx (int): x-axis length.
        ly (int): x-axis length.
        area (float): Total area.
        xval (numpy.array): Indicates where there is mass.
       
    Returns:
        A tuple of numpy.array with function and derivative values.         
    """
    fval[ind, 0] = total_area(lx, ly, area, xval)
    fval[ind, 0] -= constr_values[ind]
    dfdx[ind, :] = 100/(lx * ly) * area.reshape(1, -1)
 
    if constr_values[ind] < 0:
        fval[ind, 0] *= -1 
        dfdx[ind, :] *= -1

    return fval, dfdx

def ratio_constr(fval, dfdx, ind, constr_values, nelx, nely, coord, connect, E, v, rho, alpha_par, beta_par, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, const_func, free_ind):
    """ Calculates the function and derivative of the R Ratio.

    Args:
        fval (numpy.array): Value of the constraint function.
        dfdx (numpy.array): Value of the constraint derivative.
        ind (int): Function index in the constr_func list.
        constr_func (list): Restriction functions applied.
        constr_values (list): Values of restriction functions applied. 
        nelx (int): Number of elements on the x-axis.
        nely (int): Number of elements on the y-axis.
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio. 
        rho (float): Density.
        alpha_par (float): Damping coefficient proportional to mass.
        beta_par (float): Damping coefficient proportional to stiffness.  
        p_par (int): Penalization power to stiffness.
        q_par (int): Penalization power to mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass. 
        xval (numpy.array): Indicates where there is mass.
        disp_vector (numpy.array): Displacement.
        dyna_stif (numpy.array): Dynamic stiffness matrix.
        stif_matrix (numpy.array): Stiffness matrix.
        mass_matrix (numpy.array): Mass matrix.
        omega1_par (float): 2 * pi * frequency
        const_func (float):
        free_ind (numpy.array): DOFs free.
       
    Returns:
        A tuple of numpy.array with function and derivative values.         
    """
    fval[ind, 0], fvirg, elastic_p, kinetic_e = R_ratio(disp_vector, stif_matrix, mass_matrix, omega1_par, const_func)
    fval[ind, 0] -= constr_values[ind]
    lam_par = lambda_parameter_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, fvirg, kinetic_e, free_ind)
    dfdx[ind, :] = derivative_R(coord, connect, E, v, rho, alpha_par, beta_par, omega1_par, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, lam_par, fvirg, elastic_p, kinetic_e).reshape(nelx*nely)

    if constr_values[ind] < 0:
        fval[ind, 0] *= -1 
        dfdx[ind, :] *= -1 

    return fval, dfdx

def apply_constr(fval, dfdx, constr_func, constr_values, nelx, nely, lx, ly, coord, connect, E, v, rho, alpha_par, beta_par, p_par, q_par, x_min_k, x_min_m, area, xval, disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, const_func, free_ind):
    """ Calculates the function and derivative of the constraint functions.

    Args:
        fval (numpy.array): Value of the constraint function.
        dfdx (numpy.array): Value of the constraint derivative.
        ind (int): Function index in the constr_func list.
        constr_func (list): Restriction functions applied.
        constr_values (list): Values of restriction functions applied. 
        nelx (int): Number of elements on the x-axis.
        nely (int): Number of elements on the y-axis.
        lx (int): x-axis length.
        ly (int): x-axis length.
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio. 
        rho (float): Density.
        alpha_par (float): Damping coefficient proportional to mass.
        beta_par (float): Damping coefficient proportional to stiffness.  
        p_par (int): Penalization power to stiffness.
        q_par (int): Penalization power to mass.
        x_min_k (float): Minimum relative densities to stiffness.
        x_min_m (float): Minimum relative densities to mass.  
        area (float): Total area.
        xval (numpy.array): Indicates where there is mass.
        disp_vector (numpy.array): Displacement.
        dyna_stif (numpy.array): Dynamic stiffness matrix.
        stif_matrix (numpy.array): Stiffness matrix.
        mass_matrix (numpy.array): Mass matrix.
        omega1_par (float): 2 * pi * frequency
        const_func (float):
        free_ind (numpy.array): DOFs free.
       
    Returns:
        A tuple of numpy.array with function and derivative values.         
    """
    for i in range(len(constr_func)):
        if constr_func[i] == 'Area':
            fval, dfdx = area_constr(fval, dfdx, i, constr_values, lx, ly, area, xval)
        
        if constr_func[i] == 'R Ratio':
            fval, dfdx = ratio_constr(fval, dfdx, i, constr_values, nelx, nely, coord, connect, E, v, rho, alpha_par, beta_par, p_par, q_par, x_min_k, x_min_m, xval, disp_vector, dyna_stif, stif_matrix, mass_matrix, omega1_par, const_func, free_ind)
        
    return fval, dfdx

def update_lists(outit, fval, f0val, list_iter, list_fvals, list_f0val, constr_func, constr_values):
    """ Add new values to list of functions to plot convergence

    Args:
        outit (int): Iteration.
        fval (numpy.array): Constraint function.
        f0val (numpy.array): Objective function.
        list_iter (list): All iteration values.
        list_fvals (list): All constraint function values.
        list_f0val (list): All objective function values.
        constr_func (list): Restriction functions applied.
        constr_values (list): Values of restriction functions applied.

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
        coord (numpy.array): Coordinates of the element.
        ind (numpy.array): Element connectivity.

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
        lx (int): X-axis length.
        ly (int): Y-axis length.
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.

    Returns:
        Total element area.
    """
    return (100/(lx * ly)) * np.sum(xval * area)

def get_neighbors_radius(nelx, nely, coord, connect, radius):
    """ Check neighboring elements that have the centroid within the predetermined radius.

    Args:
        nelx (int): Number of elements on the x axis.
        nely (int): Number of elements on the x axis
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        radius (float): Radius to get elements in the vicinity of each element.

    Returns:
        A tuple of sparse matrices with neighbors and H (radius - distance).
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
        aux = len(distance[mask])
        aux1 = np.repeat(el, aux).tolist()
        aux2 = np.arange(0, aux)
        data.extend(distance[mask])
        ind_rows.extend(aux1)
        ind_cols.extend(aux2)
        #
        if aux > cols:
            cols = aux

    H = csc_matrix((data, (ind_rows, ind_cols)), shape=(nelx*nely, cols))
    neighbors = csc_matrix((neighbors, (ind_rows, ind_cols)), shape=(nelx*nely, cols))

    return neighbors, H

def calc_xnew(H, neighbors, xval):
    """ Recalculate xval.

    Args:
        H (csc_matrix): Radius subtracted from the distance between the element and the neighbors.
        neighbors (csc_matrix): Neighbors of each element.
        xval (numpy.array): Indicates where there is mass.
        
    Returns:
        New xval values.
    """
    a = 1/np.sum(H, axis=1)
    b = np.sum(H.multiply(xval[neighbors.toarray().flatten()].reshape(H.shape)), axis=1)

    xe = np.multiply(a, b)

    return np.asarray(xe)

def get_natural_freq(range_freq, delta, disp_vector):
    """ Get the frequency with the maximum displacement.

    Args:
        range_freq (list): The initial and final frequency.
        delta (int): The step between each calculation of the displacement.
        disp_vector (numpy.array): Displacement.
    
    Returns:
        A tuple with the minimum and maximum frequency.
    """
    x = np.arange(range_freq[0], range_freq[1] + 1, delta)
    y = disp_vector
    return  x[np.where(y == y.max())]

def set_initxval(constr_func, constr_values):
    """ Calculate the initial value of xval.

    Args:
        constr_func (list): Restriction functions applied.
        constr_values (list): Values of restriction functions applied.
    
    Returns:
        A tuple with the minimum and maximum frequency.
    """
    idx = [i for i, e in enumerate(constr_func) if e == constr_func[0]]
    if len(idx) > 1:
        initial_xval = ((abs(constr_values[idx[0]]) + abs(constr_values[idx[1]]))/2)/100
    else:
        initial_xval = abs(constr_values[idx[0]])/100

    return initial_xval

def get_first_freq(nelx, nely, lx, ly, func_name, force_matrix, restri_matrix=None, freq_rsp=[0, 400, 5], constr_func=['Area'], constr_values=[50], const_func=100, modes=None, rho=7860, E=210e9, v=0.3, x_min_k=1e-8, x_min_m=1e-12, alpha=0, beta=5e-6, eta=0, p_par=3, q_par=1, save=True):
        
    coord, connect, ind_rows, ind_cols = fc.regularmeshQ4(lx, ly, nelx, nely)
    initial_xval = set_initxval(constr_func, constr_values)
    xval = initial_xval * np.ones((nelx * nely, 1))
    load_vector = fc.get_load_vector(nelx, nely, force_matrix)

    restricted_dofs = fc.get_dofs(restri_matrix)
    free_ind = fc.remove_dofs(nelx, nely, restricted_dofs)
    ngl = 2 * ((nelx + 1) * (nely + 1))

    disp_vector = fc_opt.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, ngl, E, v, rho, alpha, beta, eta, 
                xval, x_min_k, x_min_m, p_par, q_par, freq_rsp[:2], freq_rsp[2], func_name, const_func, modes,
                load_vector, unrestricted_ind=free_ind)

    freq = get_natural_freq(freq_rsp[:2], freq_rsp[2], disp_vector)
    print(freq)
    plt_opt.freqresponse(freq_rsp[:2], freq_rsp[2], disp_vector, func_name, save)
    
    return freq