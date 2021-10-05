import functions_2d as fc

import cmath
import numpy as np
from scipy.sparse.linalg import spsolve

def lambda_local_ep(ngl, ind_passive, passive_el, disp_vector, dyna_stif, coord, connect, E, v, rho):
    """ Calculates the lambda parameter of the local elastic potential energy function.

    Args:
        ngl (:obj:`int`): Degrees of freedom.
        ind_passive (:obj:`numpy.array`): Index of passive elements.
        passive_el (:obj:`numpy.array`): Passive element nodes.
        disp_vector (:obj:`numpy.array`): Displacement vector.
        dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix.
        omega_par (:obj:`float`): 2 * pi * frequency.
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.

    Returns:
        Lambda parameter solution.
    """
    aux1 = np.zeros(ngl, dtype=complex)
    fadj = 0
    for i, el in enumerate(passive_el):
        Ke, _ = fc.matricesQ4(el, coord, connect, E, v, rho)
        aux1[ind_passive[i]] = Ke@disp_vector[ind_passive[i]].conjugate()
        fadj += aux1
        aux1[:] = 0 
    fadj *= -1/2
    lam = spsolve(dyna_stif, fadj)
    return lam

def lambda_local_ki(ngl, ind_passive, passive_el, disp_vector, dyna_stif, omega_par, coord, connect, E, v, rho):
    """ Calculates the lambda parameter of the local kinetic energy function.

    Args:
        ngl (:obj:`int`): Degrees of freedom.
        ind_passive (:obj:`numpy.array`): Index of passive elements.
        passive_el (:obj:`numpy.array`): Passive element nodes.
        disp_vector (:obj:`numpy.array`): Displacement vector.
        dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix.
        omega_par (:obj:`float`): 2 * pi * frequency.
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.

    Returns:
        Lambda parameter solution.
    """
    aux = np.zeros(ngl, dtype=complex)
    fadj = 0
    for i, el in enumerate(passive_el):
        _, Me = fc.matricesQ4(el, coord, connect, E, v, rho)
        aux[ind_passive[i]] = Me@disp_vector[ind_passive[i]].conjugate()
        fadj += aux
        aux[:] = 0
    fadj *= -(omega_par**2)/2
    lam = spsolve(dyna_stif, fadj)
    return lam

def lambda_compliance(disp_vector, load_vector, function):
    """ Calculates the lambda parameter of the compliance function.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement vector.
        load_vector (:obj:`numpy.array`): Force vector.
        function (:obj:`float`): Function value.

    Returns:
        Lambda parameter solution.
    """
    lam = (disp_vector.conjugate()@load_vector)/function
    return lam

def lambda_ep(disp_vector, stif_matrix, dyna_stif, free_ind):
    """ Calculates the lambda solution of the elastic potential energy function.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement vector.
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix. 
        free_ind (:obj:`numpy.array`): Free dofs.

    Returns:
        Lambda parameter solution.
    """
    lam = np.zeros(stif_matrix.shape[0], dtype=complex)
    aux = -(1/2) * (stif_matrix[free_ind, :][:, free_ind]@disp_vector[free_ind].conjugate())
    lam[free_ind] = spsolve(dyna_stif[free_ind, :][:, free_ind], aux)
    return lam

def lambda_ek(disp_vector, mass_matrix, dyna_stif, omega_par, free_ind):
    """ Calculates the lambda solution of the kinetic energy function.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement vector.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        dyna_stif (array): Stifness matrix.
        omega_par (:obj:`float`): 2 * pi * frequency.
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

def lambda_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega_par, fvirg, kinetic_e, free_ind):
    """ Calculates the lambda solution of the strain-to-kinetic function.

    Args:
        disp_vector (:obj:`numpy.array`): Displacement vector.
        dyna_stif (array): Stifness matrix. 
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        omega_par (:obj:`float`): 2 * pi * frequency.
        fvirg (:obj:`float`): Strain-to-kinetic function.
        kinetic_e (:obj:`float`):: Kinetic energy.
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

# Or use @nb.njit
# @nb.jit(nopython=True)
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
        disp_vector (:obj:`numpy.array`): Displacement vector.
        lam (:obj:`float`): Lambda parameter.

    Returns:
        Derivative of the compliance function.
    """
    deriv_f = np.empty((len(connect), 1))
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

def derivative_input_power(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector):
    """ Calculates the derivative of the input power function.
    
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
        disp_vector (:obj:`numpy.array`): Displacement vector.
        
    Returns:
        Derivative of the input power function.
    """
    deriv_f = np.empty((len(connect), 1))
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
    return deriv_f

def derivative_ep(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam):
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
        fvirg (:obj:`float`): Elastic potential energy function.

    Returns:
        Derivative elastic potential energy function.
    """
    deriv_ep = np.empty((len(connect), 1), dtype=complex)
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
    return deriv_ep

def derivative_ek(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam):
    """ Calculates the derivative of the kinetic energy function.

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
        disp_vector (:obj:`numpy.array`): Displacement vector.
        lam (:obj:`float`): Lambda parameter.

    Returns:
        Derivative of the kinetic energy function.
    """
    deriv_ek = np.empty((len(connect), 1), dtype=complex)
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
    return deriv_ek

def derivative_R(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam, fvirg, kinetic_e):
    """ Calculates the derivative of the strain-to-kinetic function.

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
        disp_vector (:obj:`numpy.array`): Displacement vector.
        lam (:obj:`float`): Lambda parameter.
        fvirg (:obj:`float`): Strain-to-kinetic function.
        kinetic_e (:obj:`float`): Kinetic energy function.

    Returns:
        Derivative of the strain-to-kinetic function function.
    """
    deriv_R = np.empty((len(connect), 1), dtype=complex)
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
        
        deriv_R[el, 0] = 1/(4*kinetic_e) * (disp_vector[ind].conjugate()@(dKe - (omega_par**2)*fvirg*dMe)@disp_vector[ind]).real + \
                        (lam[ind]@dKed@disp_vector[ind]).real
    return deriv_R

def derivative_local_ep(passive_el, lam, ind_dofs, xval, disp_vector, connect, coord, E, v, rho, x_min_k, x_min_m, omega_par, alpha, beta, p_par, q_par):
    """ Calculates the derivative of the local elastic potential energy function.

    Args:
        passive_el (:obj:`numpy.array`): Passive element nodes.
        lam (:obj:`float`): Lambda parameter.
        ind_dofs (:obj:`numpy.array`, optional): TODO
        xval (:obj:`numpy.array`): Indicates where there is mass.
        disp_vector (:obj:`numpy.array`): Displacement vector.
        connect (:obj:`numpy.array`): Element connectivity.
        coord (:obj:`numpy.array`): Coordinates of the element.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio. 
        rho (:obj:`float`): Density.
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency.
        alpha (:obj:`float`): Damping coefficient proportional to mass.
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        p_par (:obj:`int`): Penalization power to stiffness.
        q_par (:obj:`int`): Penalization power to mass.       
        
    Returns:
        Derivative of the local elastic potential energy function.
    """
    deriv_f = np.empty((len(connect), 1), dtype=complex)

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
        
        if el in passive_el:
            deriv_f[el, 0] = (1/4) *  ((disp_vector[ind].reshape(1, -1).conjugate() @ dKe @ disp_vector[ind]) + (lam[ind].reshape(1, -1) @ dKed @ disp_vector[ind]).real)[0]
        else:
            deriv_f[el, 0] = ((lam[ind].reshape(1, -1) @ dKed @ disp_vector[ind]).real)[0]
    return deriv_f

def derivative_local_ki(coord, connect, E, v, rho, alpha, beta, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam, ind_dofs, passive_el):
    """ Calculates the derivative of the local kinetic energy function.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio. 
        rho (:obj:`float`): Density.
        alpha (:obj:`float`): Damping coefficient proportional to mass.
        beta (:obj:`float`): Damping coefficient proportional to stiffness.
        omega_par (:obj:`float`): 2 * pi * frequency.
        p_par (:obj:`int`): Penalization power to stiffness.
        q_par (:obj:`int`): Penalization power to mass. 
        x_min_m (:obj:`float`): Minimum relative densities to mass. 
        x_min_k (:obj:`float`): Minimum relative densities to stiffness.
        xval (:obj:`numpy.array`): Indicates where there is mass.
        disp_vector (:obj:`numpy.array`): Displacement vector.        
        lam (:obj:`float`): Lambda parameter.
        ind_dofs (:obj:`numpy.array`, optional): TODO
        passive_el (:obj:`numpy.array`): Passive element nodes.
        
    Returns:
        Derivative of the local input power energy function.
    """ 
    deriv_ek = np.empty((len(connect), 1), dtype=complex)
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

        if el in passive_el:
            deriv_ek[el, 0] = ((omega_par**2)/4) * disp_vector[ind].conj().reshape(1, -1)@dMe@disp_vector[ind] + (lam[ind].T@dKed@disp_vector[ind]).real   
        else:
            deriv_ek[el, 0] = (lam[ind]@dKed@disp_vector[ind]).real
    return deriv_ek

def derivative_local_R(df_ep, df_ki, fvirg):
    """ Calculates the derivative of the local strain-to-kinetic function.

    Args:
        df_ep (:obj:`numpy.array`): Elastic potential energy derivative.
        df_ki (:obj:`numpy.array`): Kinetic energy derivative.
        fvirg (:obj:`float`): Local strain-to-kinetic function.

    Returns:
        Derivative of the local strain-to-kinetic function function.
    """
    #fvirg = (ep,ki)
    return df_ep * (1/fvirg[1]) - (fvirg[0]/fvirg[1]**2)*df_ki

def derivatives_objective(func_name, fvirg, disp_vector, coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, load_vector=None, mass_matrix=None, stif_matrix=None, dyna_stif=None, free_ind=None, ind_dofs=None, ngl=None, ind_passive=None, passive_el=None):
    """ Calculates the derivative of the specified function.

    Args:
        func_name (:obj:`str`): Objective function used.
        fvirg (:obj:`float`): Non-logarithm function value.
        disp_vector (:obj:`numpy.array`): Displacement vector.
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
        xval (:obj:`numpy.array`): Indicates where there is mass.
        load_vector (:obj:`numpy.array`, optional): Force vector.
        mass_matrix (:obj:`numpy.array`, optional): Mass matrix.
        stif_matrix (:obj:`numpy.array`, optional): Stiffness matrix.
        dyna_stif (:obj:`numpy.array`, optional): Dynamic stiffness matrix.
        
        free_ind (:obj:`numpy.array`, optional): Free dofs.
        ind_dofs (:obj:`numpy.array`, optional): Defaults to None.
        ngl (:obj:`int`): Degrees of freedom.
        ind_passive (:obj:`numpy.array`, optional): Index of passive elements.
        passive_el (:obj:`numpy.array`, optional): Passive element nodes.

    Returns:
        Derivative of the specified function.
    """
    if func_name == "compliance":
        lam_par = lambda_compliance(disp_vector, load_vector, fvirg)
        df0dx   = derivative_compliance(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam_par)
    
    elif func_name == "elastic_potential_energy":
        lam_par = lambda_ep(disp_vector, stif_matrix, dyna_stif, free_ind)
        df0dx   = derivative_ep(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam_par)
        #Log Scale
        df0dx[:, 0] = 10.0*df0dx[:, 0]*np.log10(np.exp(1))/fvirg

    elif func_name == "input_power":
        df0dx = derivative_input_power(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector)
        #Log Scale
        df0dx[:, 0] = 10.0*df0dx[:, 0]*np.log10(np.exp(1))/fvirg 

    elif func_name == "kinetic_energy":
        lam_par = lambda_ek(disp_vector, mass_matrix, dyna_stif, omega_par, free_ind)
        df0dx   = derivative_ek(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam_par)
        #Log Scale
        df0dx[:, 0] = 10.0*df0dx[:, 0]*np.log10(np.exp(1))/fvirg

    elif func_name == "r_ratio":
        if omega_par == 0:
            omega_par = 1e-12
        kinetic_e = ((1/4) * omega_par**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
        lam_par   = lambda_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega_par, fvirg, kinetic_e, free_ind)
        df0dx     = derivative_R(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam_par, fvirg, kinetic_e)
        #Log Scale
        df0dx[:, 0] = 10.0*df0dx[:, 0]*np.log10(np.exp(1))/fvirg

    elif func_name == "local_ep":
        lam_par = lambda_local_ep(ngl, ind_passive, passive_el, disp_vector, dyna_stif, coord, connect, E, v, rho)
        df0dx = derivative_local_ep(passive_el, lam_par, ind_dofs, xval, disp_vector, connect, coord, E, v, rho, x_min_k, x_min_m, omega_par, alpha_par, beta_par, p_par, q_par)
        #Log Scale
        df0dx[:, 0] = 10*df0dx[:, 0] * np.log10(np.exp(1))/fvirg

    elif func_name == "local_ki":
        lam_par = lambda_local_ki(ngl, ind_passive, passive_el, disp_vector, dyna_stif, omega_par, coord, connect, E, v, rho)
        df0dx = derivative_local_ki(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam_par, ind_dofs, passive_el)
        #Log Scale
        df0dx[:, 0] = 10.0*df0dx[:, 0]*np.log10(np.exp(1))/fvirg

    elif func_name == "local_r":        
        lam_par = lambda_local_ep(ngl, ind_passive, passive_el, disp_vector, dyna_stif, coord, connect, E, v, rho)
        df_ep = derivative_local_ep(passive_el, lam_par, ind_dofs, xval, disp_vector, connect, coord, E, v, rho, x_min_k, x_min_m, omega_par, alpha_par, beta_par, p_par, q_par)

        lam_par = lambda_local_ki(ngl, ind_passive, passive_el, disp_vector, dyna_stif, omega_par, coord, connect, E, v, rho)
        df_ki = derivative_local_ki(coord, connect, E, v, rho, alpha_par, beta_par, omega_par, p_par, q_par, x_min_m, x_min_k, xval, disp_vector, lam_par, ind_dofs, passive_el)
        
        df0dx = derivative_local_R(df_ep, df_ki, fvirg)
        #Log Scale
        df0dx[:, 0] = 10.0*df0dx[:, 0]*np.log10(np.exp(1))/(fvirg[0]/fvirg[1])
    return df0dx.real