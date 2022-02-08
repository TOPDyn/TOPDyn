import numpy as np
from numpy.lib.function_base import gradient
from fem_opt import FemOpt
from objective_functions import ObjectiveFunctions
from derivatives import Derivatives

class Constraint:
    """ Constraint functions """ 
    def __init__(self, constr_func, constr_values, nelx, nely, lx, ly, coord, connect, E, v, rho, ind_rows, ind_cols, alpha, beta, \
                eta, x_min_k, x_min_m, p_par, q_par, free_ind, load_vector, modes, const_func, ind_passive, passive_el, gradients=True) -> None:
        """
        Args:
            constr_func (:obj:`list`)  : Constraint functions applied.
            constr_values (:obj:`list`): Values of constraint functions applied.
            nelx (:obj:`int`): Number of elements on the x-axis. 
            nely (:obj:`int`): Number of elements on the y-axis. 
            lx (:obj:`float`): X-axis length.
            ly (:obj:`float`): Y-axis length.
            coord (:obj:`numpy.array`): mesh coordinates. 
            connect (:obj:`numpy.array`): Connectivity of elements
            E (:obj:`float`, optional): Elastic modulus.
            v (:obj:`float`, optional): Poisson's ratio.
            rho (:obj:`float`, optional): Density. 
            ind_rows (:obj:`numpy.array`): Node indexes to make the assembly.
            ind_cols (:obj:`numpy.array`): Node indexes to make the assembly.
            alpha (:obj:`float`, optional): Damping coefficient proportional to mass.
            beta (:obj:`float`, optional): Damping coefficient proportional to stiffness.
            eta (:obj:`float`, optional): Damping coefficient.
            x_min_k (:obj:`float`, optional): Minimum relative densities to stiffness.
            x_min_m (:obj:`float`, optional): Minimum relative densities to mass.
            p_par (:obj:`int`, optional): Penalization power to stiffness.
            q_par (:obj:`int`, optional): Penalization power to mass.
            free_ind (:obj:`numpy.array`): Free dofs.  
            load_vector (:obj:`numpy.array`): Load.
            modes (:obj:`int`, optional): The number of eigenvalues and eigenvectors desired. 
            const_func,  TODO
            ind_passive (:obj:`numpy.array`): Index of passive elements.
            passive_el (:obj:`numpy.array`): Passive element nodes.
            gradients (:obj:`bool`): If True calculates the derivatives.Defaults to True.
        """
        self.fem_opt = FemOpt(coord, connect, E, v, rho, nelx, nely, ind_rows, ind_cols, \
                        x_min_k, x_min_m, p_par, q_par, free_ind, load_vector, modes)

        self.obj_func = ObjectiveFunctions(load_vector, const_func, ind_passive, passive_el, coord, connect, E, v, rho)

        if gradients:
            self.dif_func = Derivatives(self.fem_opt.ngl, free_ind, ind_passive, passive_el, coord, connect, E, v, rho, alpha, beta, \
                                        p_par, q_par, x_min_k, x_min_m, load_vector)
        self.alpha = alpha
        self.beta = beta
        self.eta = eta

        self.coord = coord
        self.connect = connect
        self.lx = lx
        self.ly = ly

        aux = ["compliance", "local_ep", "local_ki", "local_r"]
        self.freq_constr_bool = any(x in aux for x in constr_func)
        if self.freq_constr_bool:
            constr_values, self.freq_constr, self.ind_freq_constr, self.ind_constr2 = self.__rewrite_constr_values(constr_func, constr_values)
            self.f_scale_constr = 100*np.ones(len(constr_values))
        else:
            self.freq_constr = None
            self.ind_constr2 = None

        self.fval = np.zeros((len(constr_func), 1))
        if gradients:
            self.dfdx = np.zeros((len(constr_func), nelx * nely))
    
        self.constr_values = np.array(constr_values)
        self.constr_func = constr_func

        self.area = self.__calc_area()
    
    def calculate(self, xval, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, omega_par, gradients=True):
        """ Calculates constraint functions and their derivatives. 
        
        Args:
            xval (:obj:`numpy.array`): Indicates where there is mass.
            disp_vector (:obj:`numpy.array`): Displacement.
            dyna_stif (:obj:`numpy.array`): Dynamic stiffness matrix.
            stif_matrix (:obj:`numpy.array`): Stiffness matrix.
            mass_matrix (:obj:`numpy.array`): Mass matrix.
            damp_matrix (:obj:`numpy.array`): Damping matrix.
            omega_par (:obj:`float`): 2 pi frequency.
            gradients (:obj:`float`): If True calculates the derivatives.
        """
        for ind in range(len(self.constr_func)):
            if self.constr_func[ind] == "area":
                aux_fval = self.calc_total_area(xval)
                if gradients:
                    aux_dfdx = 100/(self.lx * self.ly) * self.area.reshape(1, -1)

            elif self.constr_func[ind] == "local_ep":
                if self.ind_constr2[ind] == 0:
                    omega_localep = 2 * np.pi * self.freq_constr[ind]
                    dyna_stif_localep = self.fem_opt.assembly_dyna_stif(omega_localep, stif_matrix, mass_matrix, damp_matrix)
                    disp_vector_localep, _ = self.fem_opt.get_disp_vector(omega_localep, stif_matrix, mass_matrix, dyna_stif_localep, self.alpha, self.beta, self.eta)
                
                aux_fval, fvirg = self.obj_func.calculate("local_ep", disp_vector_localep)
                                                            
                if gradients:
                    aux_dfdx = self.dif_func.calculate("local_ep", fvirg, disp_vector_localep, omega_localep,  xval, dyna_stif=dyna_stif_localep)                     

            elif self.constr_func[ind] == "local_ki":
                if self.ind_constr2[ind] == 0:
                    omega_localki = 2 * np.pi *  self.freq_constr[ind]
                    dyna_stif_localki = self.fem_opt.assembly_dyna_stif(omega_localki, stif_matrix, mass_matrix, damp_matrix)
                    disp_vector_localki, _ = self.fem_opt.get_disp_vector(omega_localki, stif_matrix, mass_matrix, dyna_stif_localki, self.alpha, self.beta, self.eta)
                
                aux_fval, fvirg = self.obj_func.calculate("local_ki", disp_vector_localki, omega_par=omega_localki)

                if gradients:
                    aux_dfdx = self.dif_func.calculate("local_ki", fvirg, disp_vector_localki,  omega_localki, xval, dyna_stif=dyna_stif_localki)

            elif self.constr_func[ind] == "local_r":
                if self.ind_constr2[ind] == 0:
                    omega_localR = 2 * np.pi *  self.freq_constr[ind]
                    dyna_stif_localR = self.fem_opt.assembly_dyna_stif(omega_localR, stif_matrix, mass_matrix, damp_matrix)
                    disp_vector_localR, _ = self.fem_opt.get_disp_vector(omega_localR, stif_matrix, mass_matrix, dyna_stif_localR, self.alpha, self.beta, self.eta)
                
                aux_fval, fvirg = self.obj_func.calculate("local_r", disp_vector_localR, omega_par=omega_localR)

                if gradients:
                    aux_dfdx = self.dif_func.calculate("local_r", fvirg, disp_vector_localR, omega_localR, xval, dyna_stif=dyna_stif_localR)
            
            elif self.constr_func[ind] == "compliance":
                if self.ind_constr2[ind] == 0:
                    omega_comp = 2 * np.pi *  self.freq_constr[ind]
                    dyna_stif_comp = self.fem_opt.assembly_dyna_stif(omega_comp, stif_matrix, mass_matrix, damp_matrix)
                    disp_vector_comp, _ = self.fem_opt.get_disp_vector(omega_comp, stif_matrix, mass_matrix, dyna_stif_comp, self.alpha, self.beta, self.eta)

                aux_fval, fvirg = self.obj_func.calculate("compliance", disp_vector_comp)
                if gradients:
                    aux_dfdx = self.dif_func.calculate("compliance", fvirg, disp_vector_comp, omega_comp, xval)

            elif self.constr_func[ind] == "r_ratio":
                aux_fval, fvirg = self.obj_func.calculate("r_ratio", disp_vector, stif_matrix, mass_matrix, omega_par)
                if gradients:                                               
                    aux_dfdx = self.dif_func.calculate("r_ratio", fvirg, disp_vector, omega_par, xval, \
                                                    mass_matrix=mass_matrix, stif_matrix=stif_matrix, dyna_stif=dyna_stif)

            if self.constr_values[ind] < 0:
                aux_fval *= -1
                if gradients:
                    aux_dfdx *= -1 
            
            self.fval[ind, 0] = aux_fval
            if gradients:
                self.dfdx[ind, :] = aux_dfdx[:, 0]
        #return fval, dfdx

    def calc_total_area(self, xval):
        """ Calculates the total element area.

        Args:
            xval (:obj:`numpy.array`): Indicates where there is mass.

        Returns:
            Total element area.
        """
        return (100/(self.lx * self.ly)) * np.sum(xval * self.area)

    def __calc_area(self):
        """ Calculates the total area.

        Returns:
            area (:obj:`numpy.array`): The element area.
        """ 
        ind = self.connect[:, 1:] - 1
        area = np.empty((len(ind), 1))
        area[:, 0] = abs(self.coord[ind[:, 0], 1] * self.coord[ind[:, 1], 2] + self.coord[ind[:, 3], 1] * self.coord[ind[:, 0], 2] + \
                self.coord[ind[:, 1], 1] * self.coord[ind[:, 3], 2] - (self.coord[ind[:, 3], 1] * self.coord[ind[:, 1], 2] + \
                self.coord[ind[:, 0], 1] * self.coord[ind[:, 3], 2] + self.coord[ind[:, 1], 1] * self.coord[ind[:, 0], 2]))

        area[:, 0] += abs(self.coord[ind[:, 0], 1] * self.coord[ind[:, 3], 2] + self.coord[ind[:, 2], 1] * self.coord[ind[:, 0], 2] + \
                self.coord[ind[:, 3], 1] * self.coord[ind[:, 2], 2] - (self.coord[ind[:, 2], 1] * self.coord[ind[:, 3], 2] + \
                self.coord[ind[:, 0], 1] * self.coord[ind[:, 2], 2] + self.coord[ind[:, 3], 1] * self.coord[ind[:, 0], 2]))
        area *= 0.5 
        return area

    def update_f_scale_constr(self):
        if self.freq_constr_bool:
            self.f_scale_constr[self.ind_freq_constr] = self.fval[self.ind_freq_constr, 0]

    def update_dfdx(self, dfdx):
        """ Update derivative vector.

        Args:
            dfdx (:obj:`numpy.array`): New derivative vector.
        """
        self.dfdx = dfdx

    def update_fval(self, fval):
        """ Update constraint function values.

        Args:
            fval (:obj:`numpy.array`): Newconstraint function values.
        """
        self.fval = fval

    def normalize_fval(self):
        """ Normalizes constraint functions. """
        if self.freq_constr_bool:
            self.fval[:, 0] = 100 * self.fval[:, 0]/self.f_scale_constr
        
        self.fval[:, 0] -= self.constr_values

    def __rewrite_constr_values(self, constr_func, constr_values):
        """ Separates the constraint value and the frequency of the constraint function.
        
        Args:
            constr_values (:obj:`list`): Values of constraint functions applied with frequency.
            constr_func (:obj:`list`)  : Constraint functions applied.
        
        Returns:
            constr_values (:obj:`numpy.array`): Values of constraint functions applied without frequency.
            freq_constr (:obj:`numpy.array`): Frequency of contraint functions.
            ind_freq_constr (:obj:`numpy.array`): Indexes of constraint functions with frequency.
        """
        freq_constr = np.empty(len(constr_func))
        freq_constr.fill(np.nan)
        ind_freq_constr = []
        ind_constr2 = np.zeros(len(constr_func))

        aux_c  = None
        aux_ep = None
        aux_ki = None
        aux_r = None

        for i, value in enumerate(constr_values):
            if type(value) is list:
                freq_constr[i] = value[1]
                constr_values[i] = value[0]
                ind_freq_constr.append(i)

                if constr_func[i] == "compliance":
                    if aux_c is not None and aux_c == value[1]:
                        ind_constr2[i] = 1
                    aux_c = value[1]
                elif constr_func[i] == "local_ep":
                    if aux_ep is not None and aux_ep == value[1]:
                        ind_constr2[i] = 1
                    aux_ep = value[1]
                elif constr_func[i] == "local_ki":
                    if aux_ki is not None and aux_ki == value[1]:
                        ind_constr2[i] = 1
                    aux_ki = value[1]
                elif constr_func[i] == "local_r":
                    if aux_r is not None and aux_r == value[1]:
                        ind_constr2[i] = 1
                    aux_r = value[1]
        return constr_values, freq_constr, ind_freq_constr, ind_constr2