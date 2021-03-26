import numpy as np
import pyqtgraph as pg
import functions_3d as fc
import plots_3d as plt

def main(nelx, nely, nelz, lx, ly, lz, load_matrix, restri_matrix=None, E=210e9, v=0.3, rho=7860, alpha=1e-1, beta=1e-7, eta=1e-8, factor=1e9, freq=180, freq_rsp=[], plot_type='deformed', complete=True, amp=1e9, node_plot=None, timing=False):
    '''
    Args:
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        nelz (int): Number of elements on the Z-axis.
        lx (float): X-axis length.
        ly (float): Y-axis length.
        lz (float): Z-axis length.
        load_matrix (numpy.array): The columns are respectively node, x direction, y direction, z direction, force value.
        restri_matrix (:obj:`numpy.array`, optional)= The columns are respectively node, x direction, y direction, z direction. Defaults to None. 
        E (:obj:`float`, optional): Elastic modulus. Defaults to 210e9.
        v (:obj:`float`, optional): Poisson's ratio. Defaults to 0.3. 
        rho (:obj:`float`, optional): Density. Defaults to 7860.
        alpha (:obj:`float`, optional): Damping coefficient proportional to mass. Defaults to 1e-1.
        beta (:obj:`float`, optional): Damping coefficient proportional to stiffness. Defaults to 1e-7. 
        eta (:obj:`float`, optional): Damping coefficient. Defaults to 1e-8.
        factor (:obj:`float`, optional): Factor to deform the mesh. Defaults to 1e9.
        freq (:obj:`int`, optional): Optimized frequency. Defaults to 180.
        freq_rsp (:obj:`list`, optional): If len is 3, a frequency response graph of the original and optimized structure is generated. Defaults to [].
            First value is the minimum frequency of the graph.
            Second value is the maximum frequency of the graph.
            Third value is the step between each calculation of the objective function. 
        plot_type (:obj:`str`, optional): The type of the plot. It can be 'deformed', 'non_and_deformed', 'colorful' and 'animation'. Defaults to 'deformed'.
        complete (:obj:`bool`, optional): If the arrows and cones will be plotted. Defaults to True.
        amp (:obj:`int`, optional): Amplitude to generate deformation animation. Defaults to 1e9.
        node_plot (:obj:`numpy.array`, optional): The node to plot the frequency response graph. Defaults to first element of load_matrix.
            The columns are respectively node, x direction, y direction, z direction.
            It is a 1x4 matrix.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.           
    '''

    coord, connect, ind_rows, ind_cols = fc.regularmeshH8(nelx, nely, nelz, lx, ly, lz)
    free_ind = None
    # Get free indexes
    if restri_matrix is not None:
        restricted_dofs = fc.get_dofs(restri_matrix)
        free_ind = fc.remove_dofs(nelx, nely, nelz, restricted_dofs)
    # Calculate force
    load_vector = fc.get_load_vector(nelx, nely, nelz, load_matrix)
    # Calculate displacement vector
    disp_vector = fc.solution3D(coord, connect, ind_rows, ind_cols, nelx, nely, nelz, E, v, rho, alpha, beta, eta, freq=freq, timing=timing, load_vector=load_vector, unrestricted_ind=free_ind)
    # Mesh with displacement disp_vector
    disp_vector = fc.change_U_shape(disp_vector.real)
    coord_U = fc.apply_U(disp_vector, coord, factor)
    faces = plt.free_faces(coord_U, connect)
    #
    mesh1 = None
    if plot_type == 'non_and_deformed':
        mesh = plt.build_mesh(coord_U[:, 1:], faces, timing=timing)
        mesh1 = plt.build_mesh(coord[:, 1:], faces)
        mesh1.wireframe().lineColor('navy').lineWidth(4)
    elif plot_type == 'colorful':
        values = np.sqrt(disp_vector[:, 0]**2 + disp_vector[:, 1]**2 + disp_vector[:, 2]**2)
        mesh = plt.build_mesh(coord_U[:,1:], faces, scalars=values, timing=timing)
    elif plot_type == 'animation':
        plt.animation(coord, connect, amp, disp_vector)
    else:
        mesh = plt.build_mesh(coord_U[:, 1:], faces, timing=timing)
    #
    if complete:
        normalized_coord = mesh.points()
        arrows = plt.build_arrows(load_matrix, normalized_coord, timing=timing)
        cones = plt.build_cones(restri_matrix, normalized_coord, timing=timing)
    else:
        arrows, cones = None, None
    #
    if len(freq_rsp) == 3:
        app = pg.mkQApp()
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        if node_plot is None:
            node_plot = load_matrix[0, :4].reshape(1, 4)
        disp_plot = fc.freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, nelz, E, v, rho, alpha, beta, eta, freq_rsp[:2], freq_rsp[2], node_plot, load_vector, unrestricted_ind=free_ind)
        win = plt.plot_freqresponse(freq_rsp[:2], freq_rsp[2], disp_plot.real)
    plt.plot_mesh(mesh, arrows, cones, complete, mesh1)
    if len(freq_rsp) == 3: 
        app.exec_()