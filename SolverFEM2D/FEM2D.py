import numpy as np
import matplotlib.pyplot as plt
import plots_FEM2D as plt_fem
import functions2D as fc

def main(nelx, nely, lx, ly, load_matrix, restri_matrix=None, E=210e9, v=0.3, rho=7860, alpha=0, beta=0, eta=0, factor=1000, freq=180, node_plot=None, freq_rsp=[], save=False, timing=False):
    """ 
    Args:
        nelx (int): Number of elements on the x-axis.
        nely (int): Number of elements on the y-axis.
        lx (int): x-axis length.
        ly (int): y-axis length.
        load_matrix (numpy.array): The columns are respectively node, x direction, y direction, force value.
        restri_matrix (:obj:`numpy.array`, optional)= The columns are respectively node, x direction, y direction. Defaults to None. 
        E (:obj:`float`, optional): Elastic modulus. Defaults to 210e9.
        v (:obj:`float`, optional): Poisson's ratio. Defaults to 0.3. 
        rho (:obj:`float`, optional): Density. Defaults to 7860.
        alpha (:obj:`float`, optional): Damping coefficient proportional to mass. Defaults to 0.
        beta (:obj:`float`, optional): Damping coefficient proportional to stiffness. Defaults to 0. 
        eta (:obj:`float`, optional): Damping coefficient. Defaults to 0.
        factor (:obj:`float`, optional): Factor to deform the mesh. Defaults to 1000.
        freq (:obj:`int`, optional): Optimized frequency. Defaults to 180.
        node_plot (:obj:`numpy.array`, optional): The node to plot the frequency response graph. Defaults to first element of load_matrix.
            The columns are respectively node, x direction, y direction.
            It is a 1x3 matrix.
        freq_rsp (:obj:`list`, optional): If len is 3, a frequency response graph of the original and optimized structure is generated. Defaults to [].
            First value is the minimum frequency of the graph.
            Second value is the maximum frequency of the graph.
            Third value is the step between each calculation of the objective function. 
        save (:obj:`bool`, optional): if True save the optimization and frequency response graphs as PNG. Defaults to False.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.
    """
    coord, connect, ind_rows, ind_cols = fc.regularmeshQ4(lx, ly, nelx, nely, timing=timing)
    # Get free indexes
    free_ind = None
    if restri_matrix is not None:
        restricted_dofs = fc.get_dofs(restri_matrix)
        free_ind = fc.remove_dofs(nelx, nely, restricted_dofs)
    # Calculate force
    load_vector = fc.get_load_vector(nelx, nely, load_matrix)
    # Calculate displacement vector
    disp_vector = fc.solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, freq=freq, timing=timing, load_vector=load_vector, unrestricted_ind=free_ind)
    # Mesh with displacement
    disp_vector = fc.change_U_shape(disp_vector.real)
    coord_U = fc.apply_U(disp_vector, coord, factor)
    collection = plt_fem.build_collection(coord_U, connect[:, 1:])
    if len(freq_rsp) == 3:
        ax1 = plt_fem.plot_collection(lx, ly, coord_U, collection, load_matrix, restri_matrix, save=save)
        freq_range = freq_rsp[:2]
        delta = freq_rsp[2]
        if node_plot is None:
            node_plot = load_matrix[0, :3].reshape(1, 3)
        vector_U = fc.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, freq_range, delta, node_plot, load_vector, unrestricted_ind=free_ind)
        ax2 = plt_fem.plot_freqresponse(freq_range, delta, vector_U.real, save=save)
        plt.show()
    else:
        ax = plt_fem.plot_collection(lx, ly, coord_U, collection, load_matrix, restri_matrix, save=save)
        plt.show()