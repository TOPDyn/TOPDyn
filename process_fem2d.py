
import os
import sys
import ast
import numpy as np

from boundary_conditions import BoundConditions
from mesh_process import Mesh
from fem_2d import Fem2D

folder_name = 'temp'
directory = os.path.join(os.path.dirname(__file__), folder_name)
# READ PARAMS
param_file = os.path.join(directory, 'param_file.txt')
file = open(param_file, "r")
contents = file.read()
param = ast.literal_eval(contents)
file.close()
# READ LOAD MATRIX
load_file = os.path.join(directory, 'param_load_matrix.txt')
load_matrix = np.loadtxt(load_file)
if len(load_matrix.shape) == 1:
    load_matrix = load_matrix.reshape(1, -1)
# READ NODE CONSTRAIN MATRIX
node_constrain_file = os.path.join(directory, 'param_node_constrain_matrix.txt')
node_constrain_matrix = np.loadtxt(node_constrain_file)
if len(node_constrain_matrix.shape) == 1:
    node_constrain_matrix = node_constrain_matrix.reshape(1, -1)

# PROCESS
mesh_2d = Mesh(False, None, param["mesh_file"], param["nelx"], param["nely"], None, param["lx"], param["ly"], None)

sys.stderr.write("Total complete: 10%\n")
sys.stderr.flush()

bc = BoundConditions(False, mesh_2d.nelx, mesh_2d.nely, None, mesh_2d.coord, load_matrix, node_constrain_matrix)

sys.stderr.write("Total complete: 20%\n")
sys.stderr.flush()

fem_2d = Fem2D(mesh_2d.coord, mesh_2d.connect, param["E"], param["v"], param["rho"], mesh_2d.nelx, mesh_2d.nely, mesh_2d.ind_rows,
            mesh_2d.ind_cols, param["freq"], param["alpha"], param["beta"], param["eta"], bc.free_ind, bc.load_vector)

sys.stderr.write("Total complete: 40%\n")
sys.stderr.flush()

fem_2d.assembly_matrices()

sys.stderr.write("Total complete: 50%\n")

disp_vector = fem_2d.harmonic_solution(param["freq"])

sys.stderr.write("Total complete: 60%\n")
sys.stderr.flush()

if param["freqrsp"]:
    if param["node_plot"] is None:
        node_plot = np.array(bc.load_matrix[0,0], 0, 1, dtype='int').reshape(1, 3)
    else:
        aux = bc.get_nodes_by_coord([param["node_plot"][:2]])
        node_plot = np.array([aux[0], param["node_plot"][2], param["node_plot"][3]], dtype='int').reshape(1, 3)

    force_ind = bc.get_dofs(node_plot)
    vector_U = fem_2d.calc_freq_rsp(param["freq_range"], force_ind)

    sys.stderr.write("Total complete: 70%\n")
    sys.stderr.flush()

    dir_file = os.path.join(directory, 'node_plot.txt')
    np.savetxt(dir_file, node_plot)

    dir_file = os.path.join(directory, 'vector_U.txt')
    np.savetxt(dir_file, vector_U)

dir_file = os.path.join(directory, 'mesh_data.txt')
# open file for writing
f = open(dir_file,"w")
# write file
param = {"lx":mesh_2d.lx, "ly": mesh_2d.ly}
f.write(str(param))
# close file
f.close()

dir_file = os.path.join(directory, 'load_matrix.txt')
np.savetxt(dir_file, bc.load_matrix)

dir_file = os.path.join(directory, 'constr_matrix.txt')
np.savetxt(dir_file, bc.constr_matrix)

dir_file = os.path.join(directory, 'coord.txt')
np.savetxt(dir_file, mesh_2d.coord)

dir_file = os.path.join(directory, 'connect.txt')
np.savetxt(dir_file, mesh_2d.connect)

dir_file = os.path.join(directory, 'disp_vector.txt')
np.savetxt(dir_file, disp_vector)

sys.stderr.write("Total complete: 80%\n")
sys.stderr.flush()

# except: 
#     sys.stdout.write("Something went wrong")


# # Deformed mesh
# plot_2d = PlotsFem2d(mesh_2d.coord, mesh_2d.connect)
# disp_vector = plot_2d.change_disp_shape(disp_vector.real)

# coord_U = plot_2d.apply_disp(disp_vector, param["factor"])
# collection = plot_2d.build_collection(coord_U)




# self.parent.ax_mesh = self.parent.canvas_mesh.figure.subplots()
# plot_2d.plot_collection(self.parent.ax_mesh, mesh_2d.lx, mesh_2d.ly, coord_U, collection, bc.load_matrix, bc.constr_matrix)

# if self.parent.stop_thread:
#     self.parent.stop_thread = False
#     break

# self.update_mesh.emit(True)

# self.update_progess.emit(70)

# # Frequency response
# if param["freqrsp"]:
#     if param["node_plot"] is None:
#         node_plot = np.array(bc.load_matrix[0,0], 0, 1, dtype='int').reshape(1, 3)
#     else:
#         aux = bc.get_nodes_by_coord([param["node_plot"][:2]])
#         node_plot = np.array([aux[0], param["node_plot"][2], param["node_plot"][3]], dtype='int').reshape(1, 3)
    
#     if self.parent.stop_thread:
#         self.parent.stop_thread = False
#         break
    
#     self.update_progess.emit(75)

#     force_ind = bc.get_dofs(node_plot)
#     vector_U = fem_2d.calc_freq_rsp(param["freqrsp"], force_ind)

#     if self.parent.stop_thread:
#         self.parent.stop_thread = False
#         break

#     self.update_progess.emit(95)

#     self.parent.ax_freq = self.parent.canvas_mesh.figure.subplots()
#     plot_2d.plot_freq_rsp(self.parent.ax_freq, node_plot, param["freqrsp"], vector_U)

#     if self.parent.stop_thread:
#         self.parent.stop_thread = False
#         break

#     self.update_freq.emit(True)

#     self.update_progess.emit(98)

# if param["save"]:
#     self.parent.canvas_mesh.figure.savefig('mesh.png')
#     if param["freqrsp"]:
#         self.parent.canvas_freq.figure.savefig('frequency_response.png')

# self.update_progess.emit(100)
# break

#         self.complete_worker.emit(True)