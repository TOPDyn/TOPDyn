import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from global_matrices import GlobalMatrices

class Fem2D():
    def __init__(self, coord, connect, E, v, rho, nelx, nely, ind_rows, ind_cols, freq, alpha, beta, eta, free_ind, load_vector):
        
        self.global_matrices = GlobalMatrices(coord, connect, E, v, rho)

        self.nelx = nelx
        self.nely = nely
        self.ngl = 2 * ((self.nelx + 1) * (self.nely + 1))
        self.ind_rows = ind_cols 
        self.ind_cols = ind_rows

        self.freq = freq

        self.alpha = alpha
        self.beta = beta
        self.eta = eta
        
        self.free_ind = free_ind
        self.load_vector = load_vector

        # Calculate displacement vector
        self.assembly_matrices()
        self.disp_vector = self.harmonic_solution(self.freq)
        
    def assembly_matrices(self):
        """ Assembly matrices.

        Args:
            coord (:obj:`numpy.array`): Coordinates of the element.
            connect (:obj:`numpy.array`): Element connectivity.
            ind_rows (:obj:`numpy.array`): Node indexes to make the assembly.
            ind_cols (:obj:`numpy.array`): Node indexes to make the assembly.
            nelx (:obj:`int`): Number of elements on the X-axis.
            nely (:obj:`int`): Number of elements on the Y-axis.
            E (:obj:`float`): Elastic modulus.
            v (:obj:`float`): Poisson's ratio.  
            rho (:obj:`float`): Density.  

        Returns:
            Stiffness and mass matrices.
        """
        data_k = np.zeros((self.nelx * self.nely, 64), dtype=float)
        data_m = np.zeros((self.nelx * self.nely, 64), dtype=float)
        
        for el in range(self.nelx * self.nely):
            Ke, Me = self.global_matrices.matrices_Q4(el)
            data_k[el,:] = Ke.flatten() 
            data_m[el,:] = Me.flatten() 
        
        data_k = data_k.flatten()
        data_m = data_m.flatten()
        self.stif_matrix = csc_matrix((data_k, (self.ind_rows, self.ind_cols)), shape=(self.ngl, self.ngl))
        self.mass_matrix = csc_matrix((data_m, (self.ind_rows, self.ind_cols)), shape=(self.ngl, self.ngl))

    def harmonic_solution(self, freq):
        """ Direct method and no damping

        Args:
            alpha (:obj:`float`): Damping coefficient proportional to mass. 
            beta (:obj:`float`): Damping coefficient proportional to stiffness.  
            eta (:obj:`float`): Damping coefficient. 
            freq (:obj:`int`): Analyzed frequency.

        Returns:
            Displacement array.
        """
        omega = 2 * np.pi * freq
        damp_matrix = self.alpha * self.mass_matrix + (self.beta + self.eta/omega) * self.stif_matrix
        
        if self.free_ind is not None:
            F = self.load_vector[self.free_ind]
            K = self.stif_matrix[self.free_ind, :][:, self.free_ind]
            M = self.mass_matrix[self.free_ind, :][:, self.free_ind]
            C = damp_matrix[self.free_ind, :][:, self.free_ind]
            #Kd = -(omega**2)*M + K
            Kd = -(omega**2) * M + 1j * omega * C + K
            U = spsolve(Kd,F)
            disp_vector = np.zeros((self.ngl), dtype=complex)
            disp_vector[self.free_ind] = U
        else:
            Kd = -(omega**2) * self.mass_matrix + 1j * omega * damp_matrix + self.stif_matrix
            disp_vector = spsolve(Kd, self.load_vector)
        return disp_vector  
        
    def calc_freq_rsp(self, freq_rsp, force_ind):
        """ Get the displacement values for a specific node.

        Args:
            coord (:obj:`numpy.array`): Coordinates of the element.
            connect (:obj:`numpy.array`): Element connectivity.
            ind_rows (:obj:`numpy.array`): Node indexes to make the assembly.
            ind_cols (:obj:`numpy.array`): Node indexes to make the assembly.
            nelx (:obj:`int`): Number of elements on the X-axis.
            nely (:obj:`int`): Number of elements on the Y-axis.
            E (:obj:`float`): Elastic modulus.
            v (:obj:`float`): Poisson's ratio.  
            rho (:obj:`float`): Density.  
            alpha (:obj:`float`): Damping coefficient proportional to mass. 
            beta (:obj:`float`): Damping coefficient proportional to stiffness.  
            eta (:obj:`float`): Damping coefficient. 
            freq_rsp (:obj:`list`): Frequency range.
                First value is the minimum frequency.
                Second value is the maximum frequency.
                Third value is the step between each calculation of the objective function. 
            node_plot (:obj:`int`): TODO ESTOU PASSANDO OUTRA COISA
            load_vector (:obj:`numpy.array`): Load.

        Returns:
            Displacement array.        
        """
        if freq_rsp[0] == 0:
            freq_rsp[0] = 1e-12
        interval = np.arange(freq_rsp[0], freq_rsp[1] + 1, freq_rsp[2])
        vector_U = np.empty((len(interval)), dtype=complex)
         
        for n in range(len(interval)):
            disp_vector = self.harmonic_solution(interval[n])
            vector_U[n] = disp_vector[force_ind]
        return vector_U