from boundary_conditions import BoundConditions
from mesh_process import Mesh
from fem_opt import FemOpt
from filters import Filters
from objective_functions import ObjectiveFunctions
from derivatives import Derivatives
from constraint import Constraint
from passive_el import PassiveElements
from xval import Xval
from save_data import SaveData
from optimization import gcmmasub, kktcheck, asymp, concheck, raaupdate
import numpy as np
import os
import sys
import ast

def map_vector(vector, in_min, in_max, out_min, out_max):
  return (vector - in_min) * (out_max - out_min) / (in_max - in_min) + out_min

folder_name = 'temp'
directory = os.path.join(os.path.dirname(__file__), folder_name)
dir_file = os.path.join(directory, 'param_file.txt')
# Read parameters
file = open(dir_file, "r")
contents = file.read()
param = ast.literal_eval(contents)
file.close()  

# Read load matrix
load_file = os.path.join(directory, 'param_load_matrix.txt')
load_matrix = np.loadtxt(load_file)
if len(load_matrix.shape) == 1:
    load_matrix = load_matrix.reshape(1, -1)

# Read node constrain matrix
node_constrain_file = os.path.join(directory, 'param_node_constrain_matrix.txt')
node_constrain_matrix = np.loadtxt(node_constrain_file)
if len(node_constrain_matrix.shape) == 1:
    node_constrain_matrix = node_constrain_matrix.reshape(1, -1) 

# Vars
func_name2 = param["func_name2"]
multiobj_bool = (func_name2 is not None) and (param["n1"] != 1)
if multiobj_bool:
    freq2 = param["freq2"]
    omega2_par = 2 * np.pi * freq2 
omega1_par = 2 * np.pi * param["freq"]

########## 2D #############
mesh_2d = Mesh(False, None, param["mesh_file"], param["nelx"], param["nely"], None, param["lx"], param["ly"], None)

bc = BoundConditions(None, mesh_2d.nelx, mesh_2d.nely, None, mesh_2d.coord, load_matrix, node_constrain_matrix)

########## CLASSES #############
fem_opt = FemOpt(mesh_2d.coord, mesh_2d.connect, param["E"], param["v"], param["rho"], mesh_2d.nelx, mesh_2d.nely, mesh_2d.ind_rows, \
                mesh_2d.ind_cols, param["x_min_k"], param["x_min_m"], param["penal_k"], param["penal_m"], bc.free_ind, bc.load_vector, param["modes"])

filter = Filters(param["fac_ratio"], param["x_min_k"], mesh_2d.lx, mesh_2d.nelx, mesh_2d.nely, mesh_2d.coord, mesh_2d.connect, param["dens_filter"])

passive_el = PassiveElements(param["passive_coord"], filter.centroids, mesh_2d.ind_dofs)

dif_func = Derivatives(fem_opt.ngl, bc.free_ind, passive_el.ind_passive, passive_el.elements, mesh_2d.coord, mesh_2d.connect, \
                        param["E"], param["v"], param["rho"], param["alpha_par"], param["beta_par"], param["penal_k"], param["penal_m"], param["x_min_k"], param["x_min_m"], bc.load_vector)

obj_func = ObjectiveFunctions(bc.load_vector, param["const_func"], passive_el.ind_passive, passive_el.elements, mesh_2d.coord, \
                            mesh_2d.connect, param["E"], param["v"], param["rho"])

constraint = Constraint(param["constr_func"], param["constr_values"], mesh_2d.nelx, mesh_2d.nely, mesh_2d.lx, mesh_2d.ly, mesh_2d.coord, mesh_2d.connect, param["E"], param["v"], param["rho"], mesh_2d.ind_rows, \
                        mesh_2d.ind_cols, param["alpha_par"], param["beta_par"], param["eta_par"], param["x_min_k"], param["x_min_m"], param["penal_k"], param["penal_m"], bc.free_ind, bc.load_vector, \
                        param["modes"], param["const_func"], passive_el.ind_passive, passive_el.elements)

xval_gcmma = Xval(mesh_2d.nelx, mesh_2d.nely, constraint.constr_func, constraint.constr_values, param["dens_filter"], filter.H, filter.neighbors)

save_data = SaveData(constraint.constr_func, param["max_iter"], multiobj_bool, param["save"])

sys.stderr.write("Total complete: 10%\n")

# Create class plot
f = open(os.path.join(directory, 'param_plot.txt'), "w")
values_write = {"lx":mesh_2d.lx, "ly":mesh_2d.ly, "nelx":mesh_2d.nelx, "nely":mesh_2d.nely, "constr_func":constraint.constr_func}
f.write(str(values_write))
f.close()

np.savetxt(os.path.join(directory, 'constr_values.txt'), constraint.constr_values, fmt='%d')
sys.stderr.write("Creating plot\n")
sys.stderr.flush()

# Beam initial settings
m = len(constraint.constr_func)
n = mesh_2d.nelx * mesh_2d.nely
eeen = np.ones((n,1))
eeem = np.ones((m,1))
zerom = np.zeros((m,1))
xold1 = xval_gcmma.xval.copy()
xold2 = xval_gcmma.xval.copy()
xmin = 0.00 * eeen.copy()
xmax = 1 * eeen

if param["passive_coord"] is not None:
    xmin = xval_gcmma.set_passive_el_xmin(xmin, passive_el.elements)
    xval_gcmma.set_passive_el_xval(passive_el.elements)

low = xmin.copy()
upp = xmax.copy()
move = 1
c = 1000 * eeem 
d = eeem.copy()
a0 = 1
a = zerom.copy()
outeriter = 0 
kkttol = 0	
outit = 0

# GCMMA parameters
epsimin = 0.0000001
raa0 = 0.01
raa = 0.01 * eeem
raa0eps = 0.000001
raaeps = 0.000001 * eeem

# Converge criteria
kktnorm = kkttol + 10
chmax = 1
chtol = 1e-4
kconv = 0

# Calculate function values and gradients of the objective and constraints functions
if outeriter == 0:   
    xval_gcmma.calc_xnew()

    stif_matrix, mass_matrix, damp_matrix = fem_opt.assembly_matrices(xval_gcmma.xnew, param["alpha_par"], param["beta_par"])
    dyna_stif = fem_opt.assembly_dyna_stif(omega1_par, stif_matrix, mass_matrix, damp_matrix)
    disp_vector, natural_freqs = fem_opt.get_disp_vector(omega1_par, stif_matrix, mass_matrix, dyna_stif, param["alpha_par"], param["beta_par"], param["eta_par"])
    
    # Constraint function      
    constraint.calculate(xval_gcmma.xnew, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, omega1_par)

    # Objective function      
    f0val, fvirg = obj_func.calculate(param["func_name"], disp_vector, stif_matrix, mass_matrix, omega1_par)
    f0_scale = f0val

    # Derivative
    df0dx = dif_func.calculate(param["func_name"], fvirg, disp_vector, omega1_par, xval_gcmma.xnew, mass_matrix, stif_matrix, dyna_stif)

    if param["save"]:
        save_data.save_xval(outit, xval_gcmma.xnew)
        save_data.update_save_f0val1(outit, f0val, fvirg)
        save_data.update_save_fval(outit, constraint.fval)

    constraint.update_f_scale_constr()
    constraint.normalize_fval()

    # Multiobjective
    if multiobj_bool:
        dyna_stif2 = fem_opt.assembly_dyna_stif(omega2_par, stif_matrix, mass_matrix, damp_matrix)
        disp_vector2, _ = fem_opt.get_disp_vector(omega2_par, stif_matrix, mass_matrix, dyna_stif2, param["alpha_par"], param["beta_par"], param["eta_par"])
        
        # Second objective function
        f0val2, fvirg2 = obj_func.calculate(func_name2, disp_vector2, stif_matrix, mass_matrix, omega2_par)
        
        # Derivative
        df0dx2 = dif_func.calculate(func_name2, fvirg2, disp_vector2, omega2_par, xval_gcmma.xnew, mass_matrix, stif_matrix, dyna_stif2)

        if param["save"]:
            save_data.update_save_f0val2(outit, f0val2, fvirg2)
        
        # Filter
        dfdx, df0dx, df0dx2 = filter.apply_filter(xval_gcmma.xnew, constraint.dfdx, df0dx, df0dx2)       
        
        # Normalize multiobjective
        f0_scale_n2  = f0val2
        f0val2 = obj_func.normalize_objective(1 - abs(param["n1"]), f0_scale_n2, f0val2) 
        df0dx2 = dif_func.normalize_derivative(1 - abs(param["n1"]), f0_scale_n2, df0dx2)
        
        f0val = obj_func.normalize_objective(param["n1"], f0_scale, f0val)
        df0dx = dif_func.normalize_derivative(param["n1"], f0_scale, df0dx) 

        # Sum of functions and derivatives
        f0val = f0val + f0val2
        df0dx = df0dx + df0dx2
    else:
        # Filter
        dfdx, df0dx, _ = filter.apply_filter(xval_gcmma.xnew, constraint.dfdx, df0dx)   
       
        # Normalize objective f
        f0val = obj_func.normalize_objective(param["n1"], f0_scale, f0val)
        df0dx = dif_func.normalize_derivative(param["n1"], f0_scale, df0dx) 

    constraint.update_dfdx(dfdx)

    # Update plot
    f = open(os.path.join(directory, 'param_plot.txt'), "w")
    values_write = {"outeriter":outeriter, "outit":outit, "f0val":f0val}
    f.write(str(values_write))
    f.close()

    np.savetxt(os.path.join(directory, 'fval.txt'), constraint.fval)
    np.savetxt(os.path.join(directory, 'xnew_plot.txt'), xval_gcmma.xnew)
    np.savetxt(os.path.join(directory, 'xval_log.txt'), xval_gcmma.xval)
    if natural_freqs is not None:
        np.savetxt(os.path.join(directory, 'natural_freqs.txt'), natural_freqs)
    sys.stderr.write("Updating data\n")
    sys.stderr.flush()

constraint_new = Constraint(param["constr_func"], param["constr_values"], mesh_2d.nelx, mesh_2d.nely, mesh_2d.lx, mesh_2d.ly, mesh_2d.coord, mesh_2d.connect, param["E"], param["v"], param["rho"], mesh_2d.ind_rows, \
                        mesh_2d.ind_cols, param["alpha_par"], param["beta_par"], param["eta_par"], param["x_min_k"], param["x_min_m"], param["penal_k"], param["penal_m"], bc.free_ind, bc.load_vector, \
                        param["modes"], param["const_func"], passive_el.ind_passive, passive_el.elements, gradients=False)

constraint_new.update_fval(constraint.fval)
sys.stderr.write("Total complete: 20%\n")
sys.stderr.flush()
complete_vec = map_vector(np.arange(0, param["max_iter"]+1), 0, param["max_iter"], 30, 90).astype(int)
while (kktnorm > kkttol) and (outit < param["max_iter"]) and (kconv < 5):
    outit += 1
    outeriter += 1
    # The parameters low, upp, raa0 and raa are calculated:
    low,upp,raa0,raa= \
        asymp(outeriter,n,xval_gcmma.xval,xold1,xold2,xmin,xmax,low,upp,raa0,raa,raa0eps,raaeps,df0dx,constraint.dfdx)

    # The MMA subproblem is solved at the point xval:
    xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp= \
        gcmmasub(m,n,iter,epsimin,xval_gcmma.xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,constraint.fval,constraint.dfdx,a0,a,c,d)
    
    # The user should now calculate function values (no gradients) of the objective- and constraint
    # functions at the point xmma ( = the optimal solution of the subproblem).
    xval_gcmma.calc_xnew(xmma)
    
    stif_matrix, mass_matrix, damp_matrix = fem_opt.assembly_matrices(xval_gcmma.xnew, param["alpha_par"], param["beta_par"])
    dyna_stif = fem_opt.assembly_dyna_stif(omega1_par, stif_matrix, mass_matrix, damp_matrix)
    disp_vector, natural_freqs = fem_opt.get_disp_vector(omega1_par, stif_matrix, mass_matrix, dyna_stif, param["alpha_par"], param["beta_par"], param["eta_par"])
    
    # Constraint function      
    constraint_new.calculate(xval_gcmma.xnew, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, omega1_par, gradients=False)
    constraint_new.update_f_scale_constr()
    constraint_new.normalize_fval()

    # Objective function      
    f0valnew, fvirg = obj_func.calculate(param["func_name"], disp_vector, stif_matrix, mass_matrix, omega1_par)
    f0valnew = obj_func.normalize_objective(param["n1"], f0_scale, f0valnew)
    
    # Multiobjective
    if multiobj_bool:
        dyna_stif2 = fem_opt.assembly_dyna_stif(omega2_par, stif_matrix, mass_matrix, damp_matrix)
        disp_vector2, _ = fem_opt.get_disp_vector(omega2_par, stif_matrix, mass_matrix, dyna_stif2, param["alpha_par"], param["beta_par"], param["eta_par"])
        
        # Second objective function
        f0val2, fvirg2 = obj_func.calculate(func_name2, disp_vector2, stif_matrix, mass_matrix, omega2_par)
        
        # Normalization
        f0val2 = obj_func.normalize_objective(1 - abs(param["n1"]), f0_scale_n2, f0val2) 
  
        # Sum of functions and derivatives
        f0valnew += f0val2  

    # It is checked if the approximations are conservative:
    conserv = concheck(m, epsimin, f0app, f0valnew, fapp, constraint_new.fval)
    # While the approximations are non-conservative (conserv=0), repeated inner iterations are made:
    innerit = 0
    if conserv == 0:
        while (conserv == 0) and (innerit <= 15):
            innerit += 1
            # New values on the parameters raa0 and raa are calculated:
            raa0,raa = raaupdate(xmma, xval_gcmma.xval,xmin,xmax,low,upp,f0valnew,constraint_new.fval,f0app,fapp,raa0, \
                raa,raa0eps,raaeps,epsimin)
            # The GCMMA subproblem is solved with these new raa0 and raa:
            xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp = gcmmasub(m,n,iter,epsimin,xval_gcmma.xval,xmin, \
                xmax,low,upp,raa0,raa,f0val,df0dx,constraint.fval,constraint.dfdx,a0,a,c,d)
            # The user should now calculate function values (no gradients) of the objective- and constraint
            # functions at the point xmma ( = the optimal solution of the subproblem).
            xval_gcmma.calc_xnew(xmma)
            
            stif_matrix, mass_matrix, damp_matrix = fem_opt.assembly_matrices(xval_gcmma.xnew, param["alpha_par"], param["beta_par"])
            dyna_stif = fem_opt.assembly_dyna_stif(omega1_par, stif_matrix, mass_matrix, damp_matrix)
            disp_vector, natural_freqs = fem_opt.get_disp_vector(omega1_par, stif_matrix, mass_matrix, dyna_stif,  param["alpha_par"], param["beta_par"], param["eta_par"])
          
            # Constraint function      
            constraint_new.calculate(xval_gcmma.xnew, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, omega1_par, gradients=False)
            constraint_new.normalize_fval()

            # Objective function      
            f0valnew, fvirg = obj_func.calculate(param["func_name"], disp_vector, stif_matrix, mass_matrix, omega1_par)
            # Normalization
            f0valnew = obj_func.normalize_objective(param["n1"], f0_scale, f0valnew)
            
            # Multiobjective
            if multiobj_bool:
                dyna_stif2 = fem_opt.assembly_dyna_stif(omega2_par, stif_matrix, mass_matrix, damp_matrix)
                disp_vector2, _ = fem_opt.get_disp_vector(omega2_par, stif_matrix, mass_matrix, dyna_stif2, param["alpha_par"], param["beta_par"], param["eta_par"])
                
                # Second objective function
                f0val2, fvirg2 = obj_func.calculate(func_name2, disp_vector2, stif_matrix, mass_matrix, omega2_par)
                
                # Normalization
                f0val2 = obj_func.normalize_objective(1 - abs(param["n1"]), f0_scale_n2, f0val2) 
                
                # Sum of functions and derivatives
                f0valnew += f0val2   
            # It is checked if the approximations have become conservative:
            conserv = concheck(m, epsimin, f0app, f0valnew, fapp, constraint_new.fval)
        
    # Some vectors are updated:
    xold2 = xold1.copy()
    xold1 = xval_gcmma.xval.copy()
    xval_gcmma.update_xval(xmma.copy())

    # Re-calculate function values and gradients of the objective and constraints functions
    xval_gcmma.calc_xnew(xmma.copy())

    stif_matrix, mass_matrix, damp_matrix = fem_opt.assembly_matrices(xval_gcmma.xnew, param["alpha_par"], param["beta_par"])
    dyna_stif = fem_opt.assembly_dyna_stif(omega1_par, stif_matrix, mass_matrix, damp_matrix)
    disp_vector, natural_freqs = fem_opt.get_disp_vector(omega1_par, stif_matrix, mass_matrix, dyna_stif, param["alpha_par"], param["beta_par"], param["eta_par"])
    
    # Constraint function      
    constraint.calculate(xval_gcmma.xnew, disp_vector, dyna_stif, stif_matrix, mass_matrix, damp_matrix, omega1_par)

    # Objective function      
    f0val, fvirg = obj_func.calculate(param["func_name"], disp_vector, stif_matrix, mass_matrix, omega1_par)

    # Derivative
    df0dx = dif_func.calculate(param["func_name"], fvirg, disp_vector, omega1_par, xval_gcmma.xnew, mass_matrix, stif_matrix, dyna_stif)

    if param["save"]:
        save_data.save_xval(outit, xval_gcmma.xnew)
        save_data.update_save_f0val1(outit, f0val, fvirg)
        save_data.update_save_fval(outit, constraint.fval)

    constraint.normalize_fval()

    # Multiobjective
    if multiobj_bool:
        dyna_stif2 = fem_opt.assembly_dyna_stif(omega2_par, stif_matrix, mass_matrix, damp_matrix)
        disp_vector2, _ = fem_opt.get_disp_vector(omega2_par, stif_matrix, mass_matrix, dyna_stif2, param["alpha_par"], param["beta_par"], param["eta_par"])
        
        # Second objective function
        f0val2, fvirg2 = obj_func.calculate(func_name2, disp_vector2, stif_matrix, mass_matrix, omega2_par)
        
        # Derivative
        df0dx2 = dif_func.calculate(func_name2, fvirg2, disp_vector2, omega2_par, xval_gcmma.xnew, mass_matrix, stif_matrix, dyna_stif2)

        if param["save"]:
            save_data.update_save_f0val2(outit, f0val2, fvirg2)
        
        # Filter
        dfdx, df0dx, df0dx2 = filter.apply_filter(xval_gcmma.xnew, constraint.dfdx, df0dx, df0dx2)       
        
        # Normalize multiobjective
        f0val2 = obj_func.normalize_objective(1 - abs(param["n1"]), f0_scale_n2, f0val2) 
        df0dx2 = dif_func.normalize_derivative(1 - abs(param["n1"]), f0_scale_n2, df0dx2)
        
        f0val = obj_func.normalize_objective(param["n1"], f0_scale, f0val)
        df0dx = dif_func.normalize_derivative(param["n1"], f0_scale, df0dx) 

        # Sum of functions and derivatives
        f0val += f0val2
        df0dx += df0dx2
    else:
        # Filter
        dfdx, df0dx, _ = filter.apply_filter(xval_gcmma.xnew, constraint.dfdx, df0dx)   
    
        # Normalize objective f
        f0val = obj_func.normalize_objective(param["n1"], f0_scale, f0val)
        df0dx = dif_func.normalize_derivative(param["n1"], f0_scale, df0dx) 

    constraint.update_dfdx(dfdx)

    # The residual vector of the KKT conditions is calculated
    _, kktnorm, _ = \
        kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,constraint.fval,constraint.dfdx,a0,a,c,d)
    
    # conv. crit.
    if outeriter > 10: 
        chmax = max(abs(xold2 - xold1))/max(xold1)
        if chmax < chtol:
            kconv = kconv + 1
    else:
        chmax = 1

    sys.stderr.write("Total complete: "+str(complete_vec[outit])+"%\n")
    sys.stderr.flush()

    # Update plot
    f = open(os.path.join(directory, 'param_plot.txt'), "w")
    values_write = {"outeriter":outeriter, "outit":outit, "f0val":f0val}
    f.write(str(values_write))
    f.close()

    np.savetxt(os.path.join(directory, 'fval.txt'), constraint.fval)
    np.savetxt(os.path.join(directory, 'xnew_plot.txt'), xval_gcmma.xnew)
    np.savetxt(os.path.join(directory, 'xval_log.txt'), xval_gcmma.xval)
    if natural_freqs is not None:
        np.savetxt(os.path.join(directory, 'natural_freqs.txt'), natural_freqs)
    sys.stderr.write("Updating data\n")
    sys.stderr.flush()

sys.stderr.write("Total complete: 90%\n")
sys.stderr.flush()
# Frequency response
fig_freq = None
f_original = None
if param["freqrsp"]:
    print("Calculating the frequency response of the objective function")
    print('initial conditions')
    f_original  = fem_opt.calc_freq_rsp(xval_gcmma.initial_xval, param["freq_range"], param["func_name"], mesh_2d.coord, mesh_2d.connect, param["E"], param["v"], param["rho"], param["const_func"], \
                                        passive_el.ind_passive, passive_el.elements, param["alpha_plot"], param["beta_plot"], param["eta_plot"])
    print('optimized conditions')
    f_optimized = fem_opt.calc_freq_rsp(xval_gcmma.xnew, param["freq_range"], param["func_name"], mesh_2d.coord, mesh_2d.connect, param["E"], param["v"], param["rho"], param["const_func"], \
                                        passive_el.ind_passive, passive_el.elements, param["alpha_plot"], param["beta_plot"], param["eta_plot"])    
    
    sys.stderr.write("Total complete: 97%\n")
    sys.stderr.flush()
    # WRITE
    np.savetxt(os.path.join(directory, 'f_original.txt'), f_original.real)
    np.savetxt(os.path.join(directory, 'f_optimized.txt'), f_optimized.real)

    sys.stderr.write("Ploting frequency response\n")
    sys.stderr.flush()

# Deformed mesh
if param["mesh_deform"]:
    dir_file = os.path.join(directory, 'load_matrix.txt')
    np.savetxt(dir_file, bc.load_matrix)

    dir_file = os.path.join(directory, 'constr_matrix.txt')
    np.savetxt(dir_file, bc.constr_matrix)

    dir_file = os.path.join(directory, 'coord.txt')
    np.savetxt(dir_file, mesh_2d.coord)

    dir_file = os.path.join(directory, 'connect.txt')
    np.savetxt(dir_file, mesh_2d.connect)

    dir_file = os.path.join(directory, 'disp_vector.txt')
    np.savetxt(dir_file, disp_vector.real)

    sys.stderr.write("Ploting deformed mesh\n")
    sys.stderr.flush()
    
if param["save"]:
    save_data.save_data(multiobj_bool, outit, f_original, f_optimized)
    sys.stderr.write("Saving figures\n")
    sys.stderr.flush()

