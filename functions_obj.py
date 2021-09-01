import functions_2d as fc

import cmath
import numpy as np

def R_ratio_local(omega_par, disp_vector, passive_el, ind_passive, coord, connect, E, v, rho, aux_R=True):
    _, ep = elastic_potential_local(disp_vector, passive_el, ind_passive, coord, connect, E, v, rho)

    _, ki = kinetic_local(omega_par, disp_vector, passive_el, ind_passive, coord, connect, E, v, rho)

    f = (ep/ki).real
    if aux_R: #freqrsp
        fvirg = (ep,ki)
    else:
        fvirg = ep/ki
    #Log Scale
    f = 100 + 10*np.log10(f)
    return f, fvirg

def elastic_potential_local(disp_vector, passive_el, ind_passive, coord, connect, E, v, rho):
    ep2 = 0
    for i, ind_el in enumerate(ind_passive):
        Ke, _ = fc.matricesQ4(passive_el[i], coord, connect, E, v, rho)
        aux = disp_vector[ind_el].reshape(1, -1).conjugate()@Ke@disp_vector[ind_el]
        ep2+=aux
    fvirg = (1/4) * ep2[0].real
    #Log Scale
    f = 100 + 10*np.log10(fvirg)
    return f, fvirg

def kinetic_local(omega_par, disp_vector, passive_el, ind_passive, coord, connect, E, v, rho):
    ki = 0
    for i, ind_el in enumerate(ind_passive):
        _, Me = fc.matricesQ4(passive_el[i], coord, connect, E, v, rho)
        aux = disp_vector[ind_el].conj().reshape(1, -1)@Me@disp_vector[ind_el]
        ki+=aux
    fvirg = ((omega_par**2)/4)*ki[0].real
    #Log Scale
    f = 100 + 10*np.log10(fvirg)
    return f, fvirg

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
    fvirg = elastic_p.real
    #Log Scale
    elastic_p = const_func + 10 * np.log10(fvirg)
    return elastic_p, fvirg

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
    return kinetic_e, fvirg

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

def objective_funcs(func_name, disp_vector, stif_matrix=None, mass_matrix=None, load_vector=None, omega_par=None, const_func=None, passive_el=None, ind_passive=None, coord=None, connect=None, E=None, v=None, rho=None, aux_R=True):
    ''' Calculate objective function.

    Args:
        func_name (:obj:`str`): Objective function used.
        disp_vector (:obj:`numpy.array`): Displacement.
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        load_vector (:obj:`numpy.array`): Force.
        omega_par (:obj:`float`): 2 * pi * frequency.
        const_func (:obj:`float`):
        passive_el (:obj:`numpy.array`): Passive element nodes.
        ind_passive (:obj:`numpy.array`): Index of passive elements.
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio. 
        rho (:obj:`float`): Density.

    Returns:
        Objective function values.
    '''
    if func_name == "Compliance":
        f0val = compliance(disp_vector, load_vector)
        fvirg = f0val
    
    elif func_name == "Elastic Potential Energy":
        f0val, fvirg = elastic_potential_energy(disp_vector, stif_matrix, const_func)
    
    elif func_name == "Input Power":
        f0val, fvirg = input_power(disp_vector, load_vector, omega_par, const_func)
                  
    elif func_name == "Kinetic Energy":
        f0val, fvirg = kinetic_energy(disp_vector, mass_matrix, omega_par, const_func)
    
    elif func_name == 'R Ratio':
        f0val, fvirg = R_ratio(disp_vector, stif_matrix, mass_matrix, omega_par, const_func)

    elif func_name == 'Local Ep':
        f0val, fvirg = elastic_potential_local(disp_vector, passive_el, ind_passive, coord, connect, E, v, rho)

    elif func_name == 'Local Ki':
        f0val, fvirg = kinetic_local(omega_par, disp_vector, passive_el, ind_passive, coord, connect, E, v, rho)
    
    elif func_name == 'Local R':
        f0val, fvirg = R_ratio_local(omega_par, disp_vector, passive_el, ind_passive, coord, connect, E, v, rho, aux_R)

    import sys
    if not isinstance(f0val, float):
        sys.exit('f0val')
    if not isinstance(fvirg, float):
        sys.exit('fvirg')
    return f0val, fvirg