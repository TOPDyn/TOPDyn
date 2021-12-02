import numpy as np
from numpy.lib.function_base import gradient
from fem_opt import FemOpt
from objective_functions import ObjectiveFunctions
from derivatives import Derivatives

class Constraint:   
    def __init__(self, constr_func, constr_values, nelx, nely, lx, ly, coord, connect, E, v, rho, ind_rows, ind_cols, alpha, beta, \
                eta, x_min_k, x_min_m, p_par, q_par, free_ind, load_vector, modes, const_func, ind_passive, passive_el, gradients=True) -> None:

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
        """ Calculates constraint functions and derivative functions. """

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
                    aux_dfdx = self.dif_func.calculate("local_r", fvirg, disp_vector_localR, omega_localR, xval, dyna_stif=dyna_stif_localR, )
            
            elif self.constr_func[ind] == "compliance":
                if self.ind_constr2[ind] == 0:
                    omega_comp = 2 * np.pi *  self.freq_constr[ind]
                    dyna_stif_comp = self.fem_opt.assembly_dyna_stif(omega_comp, stif_matrix, mass_matrix, damp_matrix)
                    disp_vector_comp, _ = self.fem_opt.get_disp_vector(omega_comp, stif_matrix, mass_matrix, dyna_stif_comp, self.alpha, self.beta, self.eta)

                aux_fval, fvirg = self.obj_func.calculate("compliance", disp_vector_comp)
                if gradients:
                    aux_dfdx = self.dif_func.calculate("compliance", fvirg, disp_vector_comp, self.coord,  omega_comp, xval)

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
            lx (:obj:`float`): X-axis length.
            ly (:obj:`float`): Y-axis length.
            nelx (:obj:`int`): Number of elements on the X-axis.
            nely (:obj:`int`): Number of elements on the Y-axis.
        Returns:
            Total element area.
        """
        return (100/(self.lx * self.ly)) * np.sum(xval * self.area)

    def __calc_area(self):
        """ Calculates the total area. TODO ANTES SE CHAMAVA CALC_A
        Args:
            coord (:obj:`numpy.array`): Coordinates of the element.
            ind (:obj:`numpy.array`): Element connectivity.
        Returns:
            The element area.
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
        self.dfdx = dfdx

    def update_fval(self, fval):
        self.fval = fval

    def normalize_fval(self):
        """ Normalizes constraint functions.
        Args:
            f_scale_constr (:obj:`float`): Value of the function.
            fval (:obj:`numpy.array`): Constraint functions.
        
        Returns:
            Normalized functions.
        """
        if self.freq_constr_bool:
            self.fval[:, 0] = 100 * self.fval[:, 0]/self.f_scale_constr
        
        self.fval[:, 0] -= self.constr_values

    def __rewrite_constr_values(self, constr_func, constr_values):
        """ Separates the constraint value and the frequency of the constrain function.
        Args:
            constr_values (:obj:`list`): Values of constraint functions applied.
            constr_func (:obj:`list`)  : Constraint functions applied.
        Returns:
            constr_values, freq_constr, ind_freq_constr
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
            if type(value) is tuple:
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