import numpy as np
from scipy.sparse.linalg import spsolve
from global_matrices import GlobalMatrices

class Derivatives:
    """ Derivative of objective functions. """
    def __init__(self, ngl, free_ind, ind_passive, passive_el, coord, connect, E, v, rho, alpha, beta, p_par, q_par, x_min_k, x_min_m, load_vector) -> None:
        """
        Args:
            ngl (:obj:`int`): Degrees of freedom.
            free_ind (:obj:`numpy.array`): Free dofs.
            ind_passive (:obj:`numpy.array`): Index of passive elements.
            passive_el (:obj:`numpy.array`): Passive element nodes.
            coord (:obj:`numpy.array`): mesh coordinates. 
            connect (:obj:`numpy.array`): Connectivity of elements.
            E (:obj:`float`, optional): Elastic modulus.
            v (:obj:`float`, optional): Poisson's ratio.
            rho (:obj:`float`, optional): Density. 
            alpha (:obj:`float`, optional): Damping coefficient proportional to mass.
            beta (:obj:`float`, optional): Damping coefficient proportional to stiffness.
            p_par (:obj:`int`, optional): Penalization power to stiffness.
            q_par (:obj:`int`, optional): Penalization power to mass.
            x_min_k (:obj:`float`, optional): Minimum relative densities to stiffness.
            x_min_m (:obj:`float`, optional): Minimum relative densities to mass.
            load_vector (:obj:`numpy.array`): Load.
        """  
        self.ngl = ngl
        self.ind_passive = ind_passive
        self.passive_el = passive_el 
        self.load_vector = load_vector
        self.free_ind = free_ind
        self.alpha = alpha 
        self.beta = beta 
        self.p_par = p_par 
        self.q_par = q_par
        self.x_min_k = x_min_k
        self.x_min_m = x_min_m
        self.connect = connect
        dofs = 2
        self.ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,2]-1, dofs*connect[:,2],
                        dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,4]-1, dofs*connect[:,4]], dtype=int)-1).T
        self.global_matrices = GlobalMatrices(coord, connect, E, v, rho)

    def normalize_derivative(self, n, f0_scale, df0dx):
        """ Normalize derivative vector.
        
        Args:
            n (:obj:`float`): Weighting factor.
            f0_scale (:obj:`float`): Reference value.
            dfdx (:obj:`numpy.array`): Derivative vector.

        Returns:
            df0dx (:obj:`numpy.array`): Normalized derivative vector.
        """
        df0dx = n * 100 * df0dx/f0_scale
        return df0dx

    def __lambda_local_ep(self, disp_vector, dyna_stif):
        """ Calculates the lambda parameter of the local elastic potential energy function.
        
        Args:
            disp_vector (:obj:`numpy.array`): Displacement vector.
            dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix.

        Returns:
            lam (:obj:`numpy.array`): Lambda parameter solution.
        """
        aux1 = np.zeros(self.ngl, dtype=complex)
        fadj = 0
        for i, el in enumerate(self.passive_el):
            Ke, _ = self.global_matrices.matrices_Q4(el)
            aux1[self.ind_passive[i]] = Ke@disp_vector[self.ind_passive[i]].conjugate()
            fadj += aux1
            aux1[:] = 0 
        fadj *= -1/2
        lam = spsolve(dyna_stif, fadj)
        return lam

    def __lambda_local_ki(self, disp_vector, dyna_stif, omega_par):
        """ Calculates the lambda parameter of the local kinetic energy function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement vector.
            dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix.
            omega_par (:obj:`float`): Angular frequency.
           
        Returns:
            lam (:obj:`numpy.array`): Lambda parameter solution.
        """
        aux = np.zeros(self.ngl, dtype=complex)
        fadj = 0
        for i, el in enumerate(self.passive_el):
            _, Me = self.global_matrices.matrices_Q4(el)
            aux[self.ind_passive[i]] = Me@disp_vector[self.ind_passive[i]].conjugate()
            fadj += aux
            aux[:] = 0
        fadj *= -(omega_par**2)/2
        lam = spsolve(dyna_stif, fadj)
        return lam

    def __lambda_compliance(self, disp_vector, fvirg):
        """ Calculates the lambda parameter of the compliance function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement vector.
            fvirg (:obj:`float`): Function value without logarithm.

        Returns:
            lam (:obj:`numpy.array`): Lambda parameter solution.
        """
        lam = (disp_vector.conjugate()@self.load_vector)/fvirg
        return lam

    def __lambda_ep(self, disp_vector, stif_matrix, dyna_stif):
        """ Calculates the lambda solution of the elastic potential energy function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement vector.
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.
            dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix. 

        Returns:
            lam (:obj:`numpy.array`): Lambda parameter solution.
        """
        lam = np.zeros(stif_matrix.shape[0], dtype=complex)
        if self.free_ind is not None:
            aux = -(1/2) * (stif_matrix[self.free_ind, :][:, self.free_ind]@disp_vector[self.free_ind].conjugate())
            lam[self.free_ind] = spsolve(dyna_stif[self.free_ind, :][:, self.free_ind], aux)
        else:
            aux = -(1/2) * (stif_matrix@disp_vector.conjugate())
            lam = spsolve(dyna_stif, aux)
        return lam

    def __lambda_ek(self, disp_vector, mass_matrix, dyna_stif, omega_par):
        """ Calculates the lambda solution of the kinetic energy function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement vector.
            mass_matrix (:obj:`numpy.array`): Mass matrix.
            dyna_stif (array): Stifness matrix.
            omega_par (:obj:`float`): Angular frequency.
         
        Returns:
            lam (:obj:`numpy.array`): Lambda parameter solution.
        """
        lam = np.zeros(mass_matrix.shape[0], dtype=complex)
        if omega_par == 0:
            omega_par = 1e-12
        if self.free_ind is not None:
            aux = - (omega_par**2) * (mass_matrix[self.free_ind, :][:, self.free_ind]@disp_vector[self.free_ind].conjugate())
            lam[self.free_ind] = spsolve(dyna_stif[self.free_ind, :][:, self.free_ind], aux)
        else:
            aux = - (omega_par**2) * (mass_matrix@disp_vector.conjugate())
            lam = spsolve(dyna_stif, aux)
        return lam

    def __lambda_R(self, disp_vector, dyna_stif, stif_matrix, mass_matrix, omega_par, fvirg, kinetic_e):
        """ Calculates the lambda solution of the strain-to-kinetic function.

        Args:
            disp_vector (:obj:`numpy.array`): Displacement vector.
            dyna_stif (array): Dynamic matrix. 
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.
            mass_matrix (:obj:`numpy.array`): Mass matrix.
            omega_par (:obj:`float`): Angular frequency.
            fvirg (:obj:`float`): Strain-to-kinetic function without logarithm.
            kinetic_e (:obj:`float`): Kinetic energy.

        Returns:
            lam (:obj:`numpy.array`): Lambda parameter solution.
        """
        lam = np.zeros(mass_matrix.shape[0], dtype=complex)
        if omega_par == 0:
            omega_par = 1e-12
        if self.free_ind is not None:
            aux = - (1/(2*kinetic_e)) * ((stif_matrix[self.free_ind, :][:, self.free_ind] - (omega_par**2)*fvirg*mass_matrix[self.free_ind, :][:, self.free_ind])@disp_vector[self.free_ind].conjugate())
            lam[self.free_ind] = spsolve(dyna_stif[self.free_ind, :][:, self.free_ind], aux)
        else:
            aux = - (1/(2*kinetic_e)) * ((stif_matrix - (omega_par**2)*fvirg*mass_matrix)@disp_vector.conjugate())
            lam = spsolve(dyna_stif, aux)
        return lam

    def __derivative_compliance(self, omega_par, xval, disp_vector, lam):
        """ calculates the derivative of the compliance function.

        Args:
            omega_par (:obj:`float`): Angular frequency.
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement vector.
            lam (:obj:`float`): Lambda parameter.

        Returns:
            deriv_f (:obj:`numpy.array`): Compliance derivative.
        """
        deriv_f = np.empty((len(self.connect), 1))
        
        for el in range(len(self.connect)):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            ind = self.ind_dofs[el, :]
            dKe = self.p_par * (xval[el]**(self.p_par - 1))*(1-self.x_min_k) * Ke
            dCe = self.alpha * Me + self.beta * dKe
            if xval[el] > 0.1:
                dMe = self.q_par * (xval[el]**(self.q_par - 1))*(1-self.x_min_m) * Me
            else:
                dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-self.x_min_m) ) * Me        
            dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe
            deriv_f[el, 0] = (-lam *(disp_vector[ind].reshape(1, 8)@dKed@disp_vector[ind].reshape(8, 1)))[0,0].real
        return deriv_f

    def __derivative_input_power(self, omega_par, xval, disp_vector):
        """ Calculates the derivative of the input power function.
        
        Args:
            omega_par (:obj:`float`): Angular frequency.
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement vector.
            
        Returns:
            deriv_f (:obj:`numpy.array`): Input power derivative.
        """
        deriv_f = np.empty((len(self.connect), 1))
    
        for el in range(len(self.connect)):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            ind = self.ind_dofs[el, :]
            dKe = self.p_par * (xval[el]**(self.p_par - 1)) * (1-self.x_min_k) * Ke
            dCe = self.alpha * Me + self.beta * dKe
            if xval[el] > 0.1:
                dMe = self.q_par * (xval[el]**(self.q_par - 1))*(1-self.x_min_m) * Me
            else:
                dMe = ((9 * 3.512e7 * xval[el]**8 - 10 * 2.081e8 * xval[el]**9) * (1 - self.x_min_m)) * Me   
            dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe
            a = 1j * (disp_vector[ind].reshape(1, 8)@dKed@disp_vector[ind].reshape(8, 1))[0,0]
            deriv_f[el, 0] = -0.5 * omega_par * a.real 
        return deriv_f

    def __derivative_ep(self, omega_par, xval, disp_vector, lam):
        """ Calculates the derivative of the elastic potential energy function.

        Args:
            omega_par (:obj:`float`): Angular frequency.
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement.
            lam (:obj:`float`): Lambda parameter.

        Returns:
            deriv_ep (:obj:`numpy.array`): Elastic potential energy derivative.
        """
        deriv_ep = np.empty((len(self.connect), 1), dtype=complex)
    
        for el in range(len(self.connect)):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            ind = self.ind_dofs[el, :]
            #dKe1 = self.p_par * (xval[el]**(self.p_par - 1))*(1-self.x_min_k) * Ke.conjugate()
            dKe = self.p_par * (xval[el]**(self.p_par - 1))*(1-self.x_min_k) * Ke
            dCe = self.alpha * Me + self.beta * dKe
            if xval[el] > 0.1:
                dMe = self.q_par * (xval[el]**(self.q_par - 1)) * (1-self.x_min_m) * Me
            else:
                dMe = ((9 * 3.512e7 * xval[el]**8 - 10 * 2.081e8 * xval[el]**9) * (1-self.x_min_m) ) * Me 
            dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe             
            deriv_ep[el, 0] = (1/4) * (disp_vector[ind].conjugate()@dKe@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        return deriv_ep

    def __derivative_ek(self, omega_par, xval, disp_vector, lam):
        """ Calculates the derivative of the kinetic energy function.

        Args:
            omega_par (:obj:`float`): Angular frequency.
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement vector.
            lam (:obj:`float`): Lambda parameter.

        Returns:
            deriv_ek (:obj:`numpy.array`): Kinetic energy derivative.
        """
        deriv_ek = np.empty((len(self.connect), 1), dtype=complex)
    
        for el in range(len(self.connect)):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            ind = self.ind_dofs[el, :]
            dKe = self.p_par * (xval[el]**(self.p_par - 1)) * (1-self.x_min_k) * Ke
            dCe = self.alpha * Me + self.beta * dKe
            if xval[el] > 0.1:
                dMe = self.q_par * (xval[el]**(self.q_par - 1)) * (1-self.x_min_m) * Me
            else:
                dMe = ((9 * 3.512e7 * xval[el]**8 - 10 * 2.081e8 * xval[el]**9) * (1-self.x_min_m) ) * Me 
            dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe             
            deriv_ek[el, 0] = ((omega_par**2)/4) * (disp_vector[ind].conjugate()@dMe@disp_vector[ind]).real + (lam[ind]@dKed@disp_vector[ind]).real
        return deriv_ek

    def __derivative_R(self, omega_par, xval, disp_vector, lam, fvirg, kinetic_e):
        """ Calculates the derivative of the strain-to-kinetic function.

        Args:
            omega_par (:obj:`float`): Angular frequency.
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement vector.
            lam (:obj:`float`): Lambda parameter.
            fvirg (:obj:`float`): Strain-to-kinetic function without logarithm.
            kinetic_e (:obj:`float`): Kinetic energy function.

        Returns:
            deriv_R (:obj:`numpy.array`): Strain-to-kinetic derivative.
        """
        deriv_R = np.empty((len(self.connect), 1), dtype=complex)
    
        for el in range(len(self.connect)):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            ind = self.ind_dofs[el, :]
            dKe = self.p_par * (xval[el]**(self.p_par - 1))*(1-self.x_min_k) * Ke
            dCe = self.alpha * Me + self.beta * dKe
            if xval[el] > 0.1:
                dMe = self.q_par * (xval[el]**(self.q_par - 1)) * (1-self.x_min_m) * Me
            else:
                dMe = ((9 * 3.512e7 * xval[el]**8 - 10 * 2.081e8 * xval[el]**9) * (1-self.x_min_m) ) * Me 
            dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe        
            
            deriv_R[el, 0] = 1/(4*kinetic_e) * (disp_vector[ind].conjugate()@(dKe - (omega_par**2)*fvirg*dMe)@disp_vector[ind]).real + \
                            (lam[ind]@dKed@disp_vector[ind]).real
        return deriv_R

    def __derivative_local_ep(self, lam, xval, disp_vector, omega_par):
        """ Calculates the derivative of the local elastic potential energy function.

        Args:
            lam (:obj:`float`): Lambda parameter.
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement vector.
            omega_par (:obj:`float`): Angular frequency. 
            
        Returns:
            deriv_f (:obj:`numpy.array`): Local elastic potential energy derivative.
        """
        deriv_f = np.empty((len(self.connect), 1), dtype=complex)

        for el in range(len(self.connect)):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            ind = self.ind_dofs[el, :]
            dKe = self.p_par * (xval[el]**(self.p_par - 1))*(1-self.x_min_k) * Ke
            dCe = self.alpha * Me + self.beta * dKe
            
            if xval[el]>0.1:
                dMe = self.q_par * (xval[el]**(self.q_par - 1))*(1-self.x_min_m) * Me
            else:
                dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-self.x_min_m) ) * Me 
            dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe    
            
            if el in self.passive_el:
                deriv_f[el, 0] = (1/4) *  ((disp_vector[ind].reshape(1, -1).conjugate() @ dKe @ disp_vector[ind]) + (lam[ind].reshape(1, -1) @ dKed @ disp_vector[ind]).real)[0]
            else:
                deriv_f[el, 0] = ((lam[ind].reshape(1, -1) @ dKed @ disp_vector[ind]).real)[0]
        return deriv_f

    def __derivative_local_ki(self, omega_par, xval, disp_vector, lam):
        """ Calculates the derivative of the local kinetic energy function.

        Args:
            omega_par (:obj:`float`): Angular frequency.
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement vector.        
            lam (:obj:`float`): Lambda parameter.
            
        Returns:
            deriv_ek (:obj:`numpy.array`): Local input power energy derivative.
        """ 
        deriv_ek = np.empty((len(self.connect), 1), dtype=complex)
        for el in range(len(self.connect)):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            ind = self.ind_dofs[el, :]
            dKe = self.p_par * (xval[el]**(self.p_par - 1))*(1-self.x_min_k) * Ke
            dCe = self.alpha * Me + self.beta * dKe
            if xval[el]>0.1:
                dMe = self.q_par * (xval[el]**(self.q_par - 1))*(1-self.x_min_m) * Me
            else:
                dMe = ((9*3.512e7*xval[el]**8 - 10*2.081e8*xval[el]**9)*(1-self.x_min_m) ) * Me 
            dKed = dKe + omega_par * 1j * dCe - (omega_par**2) * dMe

            if el in self.passive_el:
                deriv_ek[el, 0] = ((omega_par**2)/4) * disp_vector[ind].conj().reshape(1, -1)@dMe@disp_vector[ind] + (lam[ind].T@dKed@disp_vector[ind]).real   
            else:
                deriv_ek[el, 0] = (lam[ind]@dKed@disp_vector[ind]).real
        return deriv_ek

    def __derivative_local_R(self, df_ep, df_ki, fvirg):
        """ Calculates the derivative of the local strain-to-kinetic function.

        Args:
            df_ep (:obj:`numpy.array`): Elastic potential energy derivative.
            df_ki (:obj:`numpy.array`): Kinetic energy derivative.
            fvirg (:obj:`float`): Local strain-to-kinetic function without logarithm.

        Returns:
            Derivative of the local strain-to-kinetic function function.
        """
        #fvirg = (ep,ki)
        return df_ep * (1/fvirg[1]) - (fvirg[0]/fvirg[1]**2)*df_ki

    def calculate(self, func_name, fvirg, disp_vector, omega_par, xval, mass_matrix=None, stif_matrix=None, dyna_stif=None):
        """ Calculte derivative of function.
        Args:
            func_name (:obj:`str`): Objective function used.
            fvirg (:obj:`float`): Function value without logarithm.
            disp_vector (:obj:`numpy.array`): Displacement vector. 
            omega_par (:obj:`float`): Angular frequency.
            xval (:obj:`numpy.array`): Indicates where there is mass. 
            mass_matrix (:obj:`numpy.array`, optional): Mass matrix. Default to None. 
            stif_matrix (:obj:`numpy.array`, optional): Stiffness matrix. Default to None.
            dyna_stif (:obj:`numpy.array`, optional): Dynamic stiffness matrix. Default to None.
            
        Returns:
            df0dx (:obj:`numpy.array`): Objective function derivative.
        """
        if func_name == "compliance":
            lam_par = self.__lambda_compliance(disp_vector, fvirg)
            df0dx   = self.__derivative_compliance(omega_par, xval, disp_vector, lam_par)
        
        elif func_name == "elastic_potential_energy":
            lam_par = self.__lambda_ep(disp_vector, stif_matrix, dyna_stif)
            df0dx   = self.__derivative_ep(omega_par, xval, disp_vector, lam_par)
            #Log Scale
            df0dx[:, 0] = 10 * df0dx[:, 0] * np.log10(np.exp(1))/fvirg

        elif func_name == "input_power":
            df0dx = self.__derivative_input_power(omega_par, xval, disp_vector)
            #Log Scale
            df0dx[:, 0] = 10 * df0dx[:, 0] * np.log10(np.exp(1))/fvirg 

        elif func_name == "kinetic_energy":
            lam_par = self.__lambda_ek(disp_vector, mass_matrix, dyna_stif, omega_par)
            df0dx   = self.__derivative_ek(omega_par, xval, disp_vector, lam_par)
            #Log Scale
            df0dx[:, 0] = 10 * df0dx[:, 0] * np.log10(np.exp(1))/fvirg

        elif func_name == "r_ratio":
            if omega_par == 0:
                omega_par = 1e-12
            kinetic_e = ((1/4) * omega_par**2 * (disp_vector.conjugate()@mass_matrix@disp_vector)).real
            lam_par   = self.__lambda_R(disp_vector, dyna_stif, stif_matrix, mass_matrix, omega_par, fvirg, kinetic_e)
            df0dx     = self.__derivative_R(omega_par, xval, disp_vector, lam_par, fvirg, kinetic_e)
            #Log Scale
            df0dx[:, 0] = 10 * df0dx[:, 0] * np.log10(np.exp(1))/fvirg

        elif func_name == "local_ep":
            lam_par = self.__lambda_local_ep(disp_vector, dyna_stif)
            df0dx = self.__derivative_local_ep(lam_par, xval, disp_vector, omega_par)
            #Log Scale
            df0dx[:, 0] = 10 * df0dx[:, 0] * np.log10(np.exp(1))/fvirg

        elif func_name == "local_ki":
            lam_par = self.__lambda_local_ki(disp_vector, dyna_stif, omega_par)
            df0dx = self.__derivative_local_ki(omega_par, xval, disp_vector, lam_par)
            #Log Scale
            df0dx[:, 0] = 10.0*df0dx[:, 0]*np.log10(np.exp(1))/fvirg

        elif func_name == "local_r":        
            lam_par = self.__lambda_local_ep(disp_vector, dyna_stif)
            df_ep = self.__derivative_local_ep(lam_par, xval, disp_vector, omega_par)

            lam_par = self.__lambda_local_ki(disp_vector, dyna_stif, omega_par)
            df_ki = self.__derivative_local_ki(omega_par, xval, disp_vector, lam_par)
            
            df0dx = self.__derivative_local_R(df_ep, df_ki, fvirg)
            #Log Scale
            df0dx[:, 0] = 10 * df0dx[:, 0] * np.log10(np.exp(1))/(fvirg[0]/fvirg[1])
        return df0dx.real