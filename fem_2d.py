import os
import numpy as np
import matplotlib.pyplot as plt
import plots_2d as plt_fem
import functions_2d as fc
from mesh_process_2d import import_mesh

def main(mesh_file, nelx, nely, lx, ly, force_matrix, restri_matrix=None, E=210e9, v=0.3, rho=7860, alpha=0, beta=0, eta=0, factor=1000, freq=180, node_plot=None, freq_rsp=[], save=False, timing=False):
    """ 
    Args:
        nelx (:obj:`int`): Number of elements on the x-axis.
        nely (:obj:`int`): Number of elements on the y-axis.
        lx (:obj:`float`): x-axis length.
        ly (:obj:`float`): y-axis length.
        force_matrix (:obj:`numpy.array`): It's a list of lists. The list can be:
            * [x_coordinate, y_coordinate, force_applied_x, force_applied_y, force_value]

            * [value_coordinate, column_to_compare, force_applied_x, force_applied_y, force_value, error_margin]

            It is possible to merge the two options. Examples:
                force_matrix = [[1, 1, 0, -1, 100]] -> Apply a negative force of modulus 100 N in the Y direction to the node at coordinate (1,1)
                force_matrix = [[0, 1, 1, 0, 200, 0.001]] -> Apply a positive force of modulus 200 N in X direction to all the nodes with x=0
                force_matrix = [[1, 1, 0, -1, 100], [0, 1, -1, 0, 200, 0.001]] -> Apply the two options above.
        restri_matrix (:obj:`numpy.array`, optional): It's a list of lists. Defaults to None. 
            * [x_coordinate, y_coordinate, constrain_disp_x, constrain_disp_y]

            * [value_coordinate, column_to_compare, constrain_disp_x, constrain_disp_y, error_margin]
        E (:obj:`float`, optional): Elastic modulus. Defaults to 210e9.
        v (:obj:`float`, optional): Poisson's ratio. Defaults to 0.3. 
        rho (:obj:`float`, optional): Density. Defaults to 7860.
        alpha (:obj:`float`, optional): Damping coefficient proportional to mass. Defaults to 0.
        beta (:obj:`float`, optional): Damping coefficient proportional to stiffness. Defaults to 0. 
        eta (:obj:`float`, optional): Damping coefficient. Defaults to 0.
        factor (:obj:`float`, optional): Factor to deform the mesh. Defaults to 1000.
        freq (:obj:`int`, optional): Optimized frequency. Defaults to 180.
        node_plot (:obj:`list`, optional): The node coordinate to plot the frequency response graph. Defaults to None.
            The elements are respectively node, x direction, y direction.
            If None the first element of the force matrix and the z direction is used.
        freq_rsp (:obj:`list`, optional): If len is 3, a frequency response graph of the original and optimized structure is generated. Defaults to [].
            First value is the minimum frequency of the graph.
            Second value is the maximum frequency of the graph.
            Third value is the step between each calculation of the objective function. 
        save (:obj:`bool`, optional): if True save the optimization and frequency response graphs as PNG. Defaults to False.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.
    """
    if mesh_file is not None:
        path = os.path.dirname(os.path.realpath(__file__)) 
        m_file = os.path.join(path, mesh_file)
        coord, connect = import_mesh(m_file)
        nelx, nely = fc.get_num_el(coord)
        lx, ly = fc.get_size_el(coord)
    else:
        coord, connect = fc.regularmeshQ4(lx, ly, nelx, nely, timing=timing)
    ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)

    # Force matrix
    force_matrix = fc.get_matrices(force_matrix, coord, True)
        
    # Get free indexes
    free_ind = None
    if restri_matrix is not None:
        restri_matrix = fc.get_matrices(restri_matrix, coord, False)
        restricted_dofs = fc.get_dofs(restri_matrix)
        free_ind = fc.remove_dofs(nelx, nely, restricted_dofs)
    
    # Calculate force
    load_vector = fc.get_load_vector(nelx, nely, force_matrix)
    
    # Calculate displacement vector
    stif_matrix, mass_matrix = fc.solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho)
    ngl = 2 * ((nelx + 1) * (nely + 1))
    disp_vector = fc.harmonic_solution(stif_matrix, mass_matrix, alpha, beta, eta, freq=freq, ngl=ngl, timing=timing, load_vector=load_vector, unrestricted_ind=free_ind)
    
    # Mesh with displacement
    disp_vector = fc.change_U_shape(disp_vector.real)
    coord_U = fc.apply_U(disp_vector, coord, factor)
    collection = plt_fem.build_collection(coord_U, connect[:, 1:])
    fig_mesh, ax1 = plt_fem.plot_collection(lx, ly, coord_U, collection, force_matrix, restri_matrix)

    if len(freq_rsp) == 3:
        if node_plot is None:
            node_plot = np.array(force_matrix[0,0], 0, 1, dtype='int').reshape(1, 3)
        else:
            aux = fc.get_nodes_by_coord(coord, [node_plot[:2]])
            node_plot = np.array([aux[0], node_plot[2], node_plot[3]], dtype='int').reshape(1, 3)
        print("Calculating the frequency response of the objective function")
        vector_U = fc.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, freq_rsp, node_plot, load_vector, unrestricted_ind=free_ind)
        fig_freq, ax2 = plt_fem.plot_freqresponse(node_plot[0,0], freq_rsp, vector_U.real)

    if save:
        folder_name = 'images'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        os.makedirs(directory, exist_ok=True)

        plt_fem.save_fig(fig_mesh, os.path.join(directory, 'mesh.png'))
        if len(freq_rsp) == 3:
            plt_fem.save_fig(fig_freq, os.path.join(directory, 'freqrsp.png'))

    print('Done!')
    plt.show()
