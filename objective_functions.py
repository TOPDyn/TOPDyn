import cmath
import numpy as np
from global_matrices import GlobalMatrices

class ObjectiveFunctions: #AQUI SERIA DO GLOBAL_MATRICES
    def __init__(self, load_vector, const_func, ind_passive, passive_el, coord, connect, E, v, rho):
        
        self.load_vector = load_vector
        self.ind_passive = ind_passive
        self.passive_el = passive_el
        
        # Optimization parameters
        self.const_func = const_func

        self.global_matrices = GlobalMatrices(coord, connect, E, v, rho)
        
    def normalize_objective(self, n, f0_scale, f0val):
        f0val = n * 100 * f0val/f0_scale
        return f0val

    def __R_ratio_local(self, omega_par, disp_vector, aux_R=True):
        """ Calculates the local strain-to-kinetic energy ratio function.

        Args:
            omega_par (:obj:`float`): 2 * pi * frequency.
            disp_vector (:obj:`numpy.array`): Displacement.
            aux_R (:obj:`bool`, optional): If True fvirg is a tuple with local kinetic energy and local elastic potential energy.
        
        Returns:
            Local strain-to-kinetic energy ratio on the logarithmic scale and the non-logarithmic strain-to-kinetic energy ratio.
        """
        _, ep = self.__elastic_potential_local(disp_vector)

        _, ki = self.__kinetic_local(omega_par, disp_vector)

        f = (ep/ki).real
        if aux_R: #freqrsp
            fvirg = (ep,ki)
        else:
            fvirg = ep/ki
        #Log Scale
        f = self.const_func + 10 * np.log10(f)
        return f, fvirg

    def __elastic_potential_local(self, disp_vector):
        """ Calculates the local elastic potential energy function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
        
        Returns:
            Local elastic potential energy on the logarithmic scale and the non-logarithmic local elastic potential energy.
        """
        ep2 = 0
        for i, ind_el in enumerate(self.ind_passive):
            Ke, _ = self.global_matrices.matrices_Q4(self.passive_el[i])
            aux = disp_vector[ind_el].reshape(1, -1).conjugate()@Ke@disp_vector[ind_el]
            ep2+=aux
        fvirg = (1/4) * ep2[0].real
        #Log Scale
        f = self.const_func + 10 * np.log10(fvirg)
        return f, fvirg

    def __kinetic_local(self, omega_par, disp_vector):
        """ Calculates the local kinetic energy function.

        Args:
            omega_par (:obj:`float`): 2 * pi * frequency.
            disp_vector (:obj:`numpy.array`): Displacement.
            
        Returns:
            Local kinetic energy on the logarithmic scale and the non-logarithmic local kinetic energy.
        """
        ki = 0
        for i, ind_el in enumerate(self.ind_passive):
            _, Me = self.global_matrices.matrices_Q4(self.passive_el[i])
            aux = disp_vector[ind_el].conj().reshape(1, -1)@Me@disp_vector[ind_el]
            ki+=aux
        fvirg = ((omega_par**2)/4) * ki[0].real
        #Log Scale
        f = self.const_func + 10 * np.log10(fvirg)
        return f, fvirg

    def __compliance(self, disp_vector):
        """ Calculates the compliance function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
            
        Returns:
            non-logarithmic compliance.
        """
        f = abs(np.dot(disp_vector, self.load_vector))
        return f

    def __input_power(self, disp_vector, omega_par):
        """ Calculates the input power function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
            omega_par (:obj:`float`): 2 * pi * frequency.
          
        Returns:
            Input power on the logarithmic scale and the non-logarithmic input power.
        """
        a = 1j * self.load_vector.conjugate()@disp_vector
        if omega_par == 0:
            omega_par = 1 #1e-12
        f = 0.5 * omega_par * a.real
        fvirg = f
        #Log Scale
        f = self.const_func + 10 * np.log10(f.real)
        return f, fvirg

    def __elastic_potential_energy(self, disp_vector, stif_matrix):
        """ Calculates the elastic potential energy function.
        
        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.

        Returns:
            Potential elastic energy on the logarithmic scale and the non-logarithmic potential elastic energy.
        """
        elastic_p = ((1/4) * (disp_vector.reshape(1, -1).conjugate()@stif_matrix@disp_vector))[0]
        fvirg = elastic_p.real
        #Log Scale
        elastic_p = self.const_func + 10 * np.log10(fvirg)
        return elastic_p, fvirg

    def __kinetic_energy(self, disp_vector, mass_matrix, omega_par):
        """ Calculates the kinetic energy function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
            mass_matrix (:obj:`numpy.array`): Mass matrix.
            omega_par (:obj:`float`): 2 * pi * frequency.
           
        Returns:
            Kinetic energy on the logarithmic scale and the non-logarithmic kinetic energy.
        """
        if omega_par == 0:
            omega_par = 1e-12
        kinetic_e = ((1/4) * omega_par**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
        fvirg = kinetic_e 
        #Log Scale
        kinetic_e  = self.const_func + 10 * np.log10(kinetic_e)
        return kinetic_e, fvirg

    def __R_ratio(self, disp_vector, stif_matrix, mass_matrix, omega_par):
        """ Calculates the strain-to-kinetic energy ratio R.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.
            mass_matrix (:obj:`numpy.array`): Mass matrix.
            omega_par (:obj:`float`): 2 * pi * frequency.
      
        Returns:
            Strain-to-kinetic energy ratio on the logarithmic scale and the non-logarithmic strain-to-kinetic energy ratio.
        """
        elastic_p = ((1/4) * (disp_vector.reshape(1, -1).conjugate()@stif_matrix@disp_vector))[0]
        if omega_par == 0:
            omega_par = 1e-12
        kinetic_e = ((1/4) * omega_par**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
        R = (elastic_p/kinetic_e)
        fvirg = R
        #Log Scale
        R  = self.const_func + 10 * np.log10(R)
        return R.real, fvirg.real

    def calculate(self, func_name, disp_vector, stif_matrix=None, mass_matrix=None, omega_par=None, aux_R=True):

        if func_name == "compliance":
            f0val = self.__compliance(disp_vector)
            fvirg = f0val
        
        elif func_name == "elastic_potential_energy":
            f0val, fvirg = self.__elastic_potential_energy(disp_vector, stif_matrix)
        
        elif func_name == "input_power":
            f0val, fvirg = self.__input_power(disp_vector, omega_par)
                    
        elif func_name == "kinetic_energy":
            f0val, fvirg = self.__kinetic_energy(disp_vector, mass_matrix, omega_par)
        
        elif func_name == "r_ratio":
            f0val, fvirg = self.__R_ratio(disp_vector, stif_matrix, mass_matrix, omega_par)

        elif func_name == "local_ep":
            f0val, fvirg = self.__elastic_potential_local(disp_vector)

        elif func_name == "local_ki":
            f0val, fvirg = self.__kinetic_local(omega_par, disp_vector)
        
        elif func_name == "local_r":
            f0val, fvirg = self.__R_ratio_local(omega_par, disp_vector, aux_R)
        return f0val, fvirg