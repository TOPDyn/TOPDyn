import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import eigh
from global_matrices import GlobalMatrices
from objective_functions import ObjectiveFunctions

class FemOpt:
    def __init__(self, coord, connect, E, v, rho, nelx, nely, ind_rows, ind_cols, x_min_k, x_min_m, p_par, q_par, free_ind, load_vector, modes):
      
        self.global_matrices = GlobalMatrices(coord, connect, E, v, rho)

        self.nelx = nelx
        self.nely = nely
        self.ngl = 2 * ((self.nelx + 1) * (self.nely + 1))
        self.ind_rows = ind_rows
        self.ind_cols = ind_cols
        self.x_min_k = x_min_k
        self.x_min_m = x_min_m
        self.p_par = p_par
        self.q_par = q_par 

        self.modes = modes

        self.free_ind = free_ind
        self.load_vector = load_vector

    def assembly_matrices(self, xval, alpha, beta):
        """ Gets values for assembling the stiffness and mass matrices.

        Args:
            xval (:obj:`numpy.array`): Indicates where there is mass.

        Returns:
            A tuple with stiffness matrix, mass matrix and time to assembly the matrices.
        """
        data_k = np.zeros((self.nelx * self.nely, 64), dtype=complex)
        data_m = np.zeros((self.nelx * self.nely, 64), dtype=complex)
        
        for el in range(self.nelx * self.nely):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            data_k[el, :] = (self.x_min_k + (xval[el]**self.p_par)*(1-self.x_min_k))* Ke.flatten()
            if xval[el]>0.1:
                data_m[el, :] = (self.x_min_m + (xval[el]**self.q_par)*(1-self.x_min_m)) * Me.flatten()
            else:
                data_m[el, :] =  (self.x_min_m + (3.512e7*xval[el]**9 - 2.081e8*xval[el]**10)*(1-self.x_min_m) ) * Me.flatten()       

        stif_matrix = csc_matrix((data_k.flatten(), (self.ind_rows, self.ind_cols)), shape=(self.ngl, self.ngl)).real
        mass_matrix = csc_matrix((data_m.flatten(), (self.ind_rows, self.ind_cols)), shape=(self.ngl, self.ngl)).real
        damp_matrix = (alpha * mass_matrix + (beta) * stif_matrix).real
        return stif_matrix, mass_matrix, damp_matrix

    def assembly_dyna_stif(self, omega_par, stif_matrix, mass_matrix, damp_matrix):
        """ Assembly the dynamic stiffness matrix.

        Args:
            omega_par (:obj:`float`): 2 pi frequency.
            mass_matrix (:obj:`numpy.array`): Mass matrix.
            damp_matrix (:obj:`numpy.array`): Damping matrix. 
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        
        Returns:
            Dynamic stiffness matrix.
        """
        return -(omega_par**2) * mass_matrix + 1j * omega_par * damp_matrix + stif_matrix

    def _harmonic_problem(self, dyna_stif):
        """ Calculates displacement.

        Args:        
            dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix. 
   
        Returns:
            A tuple with displacement vector and time.
        """
        if self.free_ind is not None:
            disp_vector = np.zeros((self.ngl), dtype=complex)
            disp_vector[self.free_ind] = spsolve(dyna_stif[self.free_ind, :][:, self.free_ind], self.load_vector[self.free_ind])
        else:
            disp_vector = spsolve(dyna_stif, self.load_vector)
        return disp_vector

    def _modal_analysis(self, mass_matrix, stif_matrix):
        """ Modal Analysis. Use eigh Scipy function.

        Args:
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.
            mass_matrix (:obj:`numpy.array`): Mass matrix.
        
        Returns:
            A tuple with natural frequencies and modes shape.
        """
        if self.free_ind is not None:
            eigen_values, eigen_vectors = eigh(stif_matrix[self.free_ind, :][:, self.free_ind].todense(), b=mass_matrix[self.free_ind, :][:, self.free_ind].todense(), subset_by_index=[0,self.modes])
        else:
            eigen_values, eigen_vectors = eigh(stif_matrix.todense(), b=mass_matrix.todense(), subset_by_index=[0,self.modes])
        
        # mod_freqmax = 3 * freqmax
        # eigen_values, eigen_vectors = eigh(stif_matrix.todense(), b=mass_matrix.todense(), subset_by_value=[-np.inf, (2*np.pi*mod_freqmax)])

        positive_real = np.absolute(np.real(eigen_values))
        natural_frequencies = np.sqrt(positive_real)/(2*np.pi)
        modal_shape = np.real(eigen_vectors)

        index_order = np.argsort(natural_frequencies)
        natural_frequencies = natural_frequencies[index_order]
        modal_shape = modal_shape[:, index_order]
        return natural_frequencies, modal_shape

    def _mode_superposition(self, natural_frequencies, modal_shape, omega_par, alpha, beta, eta, stif_matrix):    
        """ Perform an harmonic analysis through superposition method and returns the response of
            all nodes due the external or internal equivalent load. It has been implemented two
            different damping models: Viscous Proportional and Hysteretic Proportional
            Entries for Viscous Proportional Model Damping: (alpha_v, beta_v)
            Entries for Hyteretic Proportional Model Damping: (alpha_h, beta_h)
        
        Args:
            #TODO: FALTANDO AQUI
            omega_par (:obj:`float`): 2 pi frequency

        Returns:
            A tuple with displacement and time to solve the problem.
        """
        alphaV, betaV, betaH = alpha, beta, eta

        #natural_frequencies, modal_shape = modal_analysis(stif_matrix[free_ind, :][:, free_ind], mass_matrix[free_ind, :][:, free_ind], modes=modes)
        if self.free_ind is not None:
            F_aux = modal_shape.T @ self.load_vector[self.free_ind] 
        else:
            F_aux = modal_shape.T @ self.load_vector
        omega_n = 2*np.pi*natural_frequencies
        F_kg = (omega_n**2)

        F_mg =  - (omega_par**2)
        F_cg = 1j*((betaH + betaV*omega_par)*(omega_n**2) + (0. + omega_par*alphaV)) 
        data = np.divide(1, (F_kg + F_mg + F_cg))
        diag = np.diag(data)
        #disp_vector = modal_shape @ (diag @ F_aux[:,i])
        rows = stif_matrix.shape[0]
        disp_vector = np.zeros((rows), dtype=complex)
        if self.free_ind is not None:
            disp_vector[self.free_ind] = modal_shape @ (diag @ F_aux)
        else:
            disp_vector = modal_shape @ (diag @ F_aux)
        return disp_vector

    def get_disp_vector(self, omega_par, stif_matrix, mass_matrix, dyna_stif, alpha, beta, eta):
        """ Calculates displacement vector.

        Args:
            omega_par (:obj:`float`): 2 pi frequency.
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.
            mass_matrix (:obj:`numpy.array`): Mass matrix. 
            dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix.
           
        Returns:
            displacement vector, natural frequencies if using modal analysis.
        """
        if self.modes is not None:
            natural_frequencies, modal_shape = self._modal_analysis(mass_matrix, stif_matrix)
            disp_vector = self._mode_superposition(natural_frequencies, modal_shape, omega_par, alpha, beta, eta, stif_matrix)
        else:  
            disp_vector = self._harmonic_problem(dyna_stif)
            natural_frequencies = None
        return disp_vector, natural_frequencies

    def calc_freq_rsp(self, xval, freq_range, func_name, coord, connect, E, v, rho, const_func, ind_passive, passive_el, alpha, beta, eta):
        
        if freq_range[0] == 0:
            freq_range[0] = 1e-12
        
        interval = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
        func_vector = np.empty((len(interval)), dtype=complex)

        l = len(interval) + 5
      

        stif_matrix, mass_matrix, damp_matrix = self.assembly_matrices(xval, alpha, beta)
        # if modes is not None:
        #     natural_frequencies, modal_shape = self._modal_analysis(mass_matrix, stif_matrix)

        obj_func = ObjectiveFunctions(self.load_vector, const_func, ind_passive, passive_el, coord, connect, E, v, rho)
        
        for n in range(len(interval)):
            omega_par = 2 * np.pi * interval[n]
            dyna_stif = self.assembly_dyna_stif(omega_par, stif_matrix, mass_matrix, damp_matrix)
            disp_vector, _ = self.get_disp_vector(omega_par, stif_matrix, mass_matrix, dyna_stif, alpha, beta, eta)

            _, func_vector[n] = obj_func.calculate(func_name, disp_vector, stif_matrix, mass_matrix, omega_par, aux_R=False)
        return func_vector 

    def printProgressBar(self, iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """ Call in a loop to create terminal progress bar.
        `Code link <https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console>`_
        
        Args:
            iteration (:obj:`int`): current iteration.
            total (:obj:`int`): total iterations.
            prefix (:obj:`str`, optional): prefix string.
            suffix (:obj:`str`, optional): suffix string.
            decimals (:obj:`int`, optional): positive number of decimals in percent complete.
            length (:obj:`int`, optional): character length of bar.
            fill (:obj:`str`, optional): bar fill character.
            printEnd (:obj:`str`, optional): end character.
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total: 
            print()