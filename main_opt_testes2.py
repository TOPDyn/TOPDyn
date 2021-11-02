import os
import numpy as np
import functions_opt as fc_opt

mma = True

mesh_file = None

nelx, nely = 10, 5
lx, ly = 1, 0.5

rho = 7860
E = 210e9
v = 0.3

x_min_m = 1e-12
x_min_k = 1e-9
alpha_par, beta_par, eta_par = 0, 1e-6, 0 #Não fala nada do eta
alpha_plot, beta_plot, eta_plot = 0, 1e-6, 0

p_par = 3
q_par = 1

const_func = 100

# Factor applied in the radius
fac_ratio = 0.021/0.1 #raio/n.el #2.2 #2.1

# If not None is used mode superposition method
modes = None

# Create matrix of loads 
load_matrix = []
y_val = np.arange(0.23, 0.27, 0.1)
for y in y_val:
    aux = {"x_coord":1, "y_coord":y, "x_direc":0, "y_direc":-1, "force":1000}
    load_matrix.append(aux)

# Create constrain nodes matrix
constr_matrix = [{"coord":0, "axis":1, "eps":0.001, "constrain_disp_x":1, "constrain_disp_y":1}]

# Weight at objective function
n1 = 0.8

# It can be "compliance", "input_power", "elastic_potential_energy", "kinetic_energy", "r_ratio", "local_ep", "local_ki" or "local_r"
func_name = "input_power" 

# Frequency optimized for func_name
freq1 = 50 

# Tuple with func_name2 and frequency optimized for func_name2. Associated with weight (1 - n1)
multiobjective = ("compliance", 0)

# Constraint - The first function in the list is used to define the initial value of xval. "compliance", "local_ep", "local_ki", "local_r" -> (constraint value, frequency)(constraint value, frequency)
constr_func = ["area", "local_ep"]
constr_values = [50, (70, 50)]

passive_coord = ((0.4, 0.6), (0.15, 0.35)) # ((x_initial, x_final), (y_initial, y_final)) or None

# Frequency response plot
freq_rsp = [1, 1000, 1]

# If False use sensitivity filter
dens_filter = True

# If True plots the convergence graph for each iteration of the optimization
each_iter = True

# Plot mesh  
mesh_deform = False 
factor = 800
save = True
timing = False

# Method iterations
max_iter = 5

fc_opt.exe_opt(mma, mesh_file, nelx, nely, lx, ly, func_name, load_matrix, constr_matrix, freq1, constr_func, constr_values, n1, multiobjective, const_func, fac_ratio, modes, rho, E, v, x_min_m, x_min_k, alpha_par, beta_par, eta_par, alpha_plot, beta_plot, eta_plot, p_par, q_par, passive_coord, freq_rsp, dens_filter, each_iter, max_iter, mesh_deform, factor, save, timing)