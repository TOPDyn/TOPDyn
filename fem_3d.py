import os
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
import functions_3d as fc
import plots_3d as plt_3d
from mesh_process_3d import import_mesh

def main(mesh_file, nelx, nely, nelz, lx, ly, lz, load_matrix, restri_matrix=None, num_el=15, E=210e9, v=0.3, rho=7860, alpha=1e-1, beta=1e-7, eta=1e-8, factor=1e9, freq=180, freq_range=[], plot_type="deformed", complete=True, amp=1e9, node_plot=None, timing=False):
    """
    Args:
        mesh_file (:obj:`iges`): iges file.
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        nelz (int): Number of elements on the Z-axis.
        lx (:obj:`float`): X-axis length.
        ly (:obj:`float`): Y-axis length.
        lz (:obj:`float`): Z-axis length.
        load_matrix (:obj:`numpy.array`): List of dictionaries.
            The dictionary can be:
                * {"coord":value_coordinate, "axis":column_to_compare, "eps":error_margin, "x_direc":force_applied_x, "y_direc":force_applied_y, "z_direc":force_applied_z, "force":force_value}
                * {"x_coord":x_coordinate, "y_coord":y_coordinate, "apply_x":force_applied_x, "apply_y":force_applied_y, "apply_z":force_applied_z, "force":force_value}
            It's possible to merge the two options. Examples:
                * load_matrix = [{"x_coord":1, "y_coord":1, "z_coord":1, "apply_x":0, "apply_y":-1, "apply_y":0, "force":100}] -> Apply a negative force of modulus 100 N in the Y direction to the node at coordinate (1,1,1).
                * load_matrix = [{"coord":0, "axis":1,  "eps":0.001, "x_direc":1, "y_direc":0,  "z_direc":0, "force":200}] -> Apply a positive force of modulus 200 N in X direction to all the nodes with x=0.
        restri_matrix (:obj:`numpy.array`, optional): List of dictionaries. Defaults to None. 
            The dictionary can be:
                * {"x_coord":x_coordinate, "y_coord":y_coordinate, "z_coord":z_coordinate, "constrain_disp_x":constrain_disp_x, "constrain_disp_y":constrain_disp_y, "constrain_disp_z":constrain_disp_z}
                * {"coord":value_coordinate, "axis":column_to_compare, "eps":error_margin, "constrain_disp_x":constrain_disp_x, "constrain_disp_y":constrain_disp_y, "constrain_disp_z":constrain_disp_z} 
            It's possible to merge the two options. Examples:
                * restri_matrix = [{"x_coord":0, "y_coord":0.25, "y_coord":0.15, "constrain_disp_x":1, "constrain_disp_y":0,  "constrain_disp_z":0}] -> Constrain the nodes in X direction at coordinate (0, 0.25,0.15).
                * restri_matrix = [{"coord":0, "axis":1, "eps":0.001, "constrain_disp_x":0, "constrain_disp_y":1, "constrain_disp_z":0}] -> Constrain the nodes in Y direction to all the nodes with x=0.
        num_el (:obj:`int`, optional): Number of elements when passing an iges file. Defaults to 15.
        E (:obj:`float`, optional): Elastic modulus. Defaults to 210e9.
        v (:obj:`float`, optional): Poisson's ratio. Defaults to 0.3. 
        rho (:obj:`float`, optional): Density. Defaults to 7860.
        alpha (:obj:`float`, optional): Damping coefficient proportional to mass. Defaults to 1e-1.
        beta (:obj:`float`, optional): Damping coefficient proportional to stiffness. Defaults to 1e-7. 
        eta (:obj:`float`, optional): Damping coefficient. Defaults to 1e-8.
        factor (:obj:`float`, optional): Factor to deform the mesh. Defaults to 1e9.
        freq (:obj:`int`, optional): Optimized frequency. Defaults to 180.
        freq_range (:obj:`list`, optional): If len is 3, a frequency response graph of the original and optimized structure is generated. Defaults to [].
            First value is the minimum frequency of the graph.
            Second value is the maximum frequency of the graph.
            Third value is the step between each calculation of the objective function. 
        plot_type (:obj:`str`, optional): The type of the plot. It can be "deformed", "non_and_deformed", "colorful" and "animation". Defaults to "deformed".
        complete (:obj:`bool`, optional): If the arrows and cones will be plotted. Defaults to True.
        amp (:obj:`int`, optional): Amplitude to generate deformation animation. Defaults to 1e9.
        node_plot (:obj:`numpy.array`, optional): The node to plot the frequency response graph. Defaults to first element of load_matrix.
            The columns are respectively node, x direction, y direction, z direction.
            It is a 1x4 matrix.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.           
    """
    if mesh_file is not None:
        path = os.path.dirname(os.path.realpath(__file__)) 
        m_file = os.path.join(path, mesh_file)
        coord, connect = import_mesh(m_file, num_el)
        lx = max(coord[:, 1])
        ly = max(coord[:, 2])
        lz = max(coord[:, 3])
        nelx = len(coord[np.logical_and(coord[:, 2] == coord[0, 2], coord[:, 3] == coord[0, 3])]) - 1
        nely = len(coord[np.logical_and(coord[:, 1] == coord[0, 1], coord[:, 3] == coord[0, 3])]) - 1
        nelz = len(coord[np.logical_and(coord[:, 1] == coord[0, 1], coord[:, 2] == coord[0, 2])]) - 1
    else:
        coord, connect = fc.regularmeshH8(nelx, nely, nelz, lx, ly, lz)
    ngl = 3 * ((nelx + 1) * (nely + 1) * (nelz + 1))
    ind_rows, ind_cols = fc.generate_ind_rows_cols(connect)

    # Load matrix
    load_matrix = fc.get_matrices(load_matrix, coord, True)
    
    free_ind = None
    # Get free indexes
    if restri_matrix is not None:
        restri_matrix = fc.get_matrices(restri_matrix, coord, False)
        restricted_dofs = fc.get_dofs(restri_matrix)
        free_ind = fc.remove_dofs(nelx, nely, nelz, restricted_dofs)

    # Calculate load vector
    load_vector = fc.get_load_vector(ngl, load_matrix)

    # Calculate displacement vector
    stif_matrix, mass_matrix = fc.stif_mass_matrices(nelx, nely, nelz, coord, connect, ind_rows, ind_cols, E, v, rho, timing=timing)
    damp_matrix = fc.get_damp_matrix(mass_matrix, stif_matrix, freq, alpha, beta, eta)
    disp_vector = fc.get_displacement(load_vector, free_ind, stif_matrix, mass_matrix, damp_matrix, freq, ngl, timing=timing)
    
    # Mesh with displacement disp_vector
    disp_vector = fc.change_U_shape(disp_vector.real)
    coord_U = fc.apply_U(disp_vector, coord, factor)
    faces = plt_3d.free_faces(coord_U, connect)
    
    mesh1 = None
    if plot_type == "non_and_deformed":
        mesh = plt_3d.build_mesh(coord_U[:, 1:], faces, timing=timing)
        mesh1 = plt_3d.build_mesh(coord[:, 1:], faces)
        mesh1.wireframe().lineColor("navy").lineWidth(4)
    elif plot_type == "colorful":
        values = np.sqrt(disp_vector[:, 0]**2 + disp_vector[:, 1]**2 + disp_vector[:, 2]**2)
        mesh = plt_3d.build_mesh(coord_U[:,1:], faces, scalars=values, timing=timing)
    elif plot_type == "animation":
        plt_3d.animation(coord, connect, amp, disp_vector)
    else:
        mesh = plt_3d.build_mesh(coord_U[:, 1:], faces, timing=timing)
    
    if complete:
        normalized_coord = mesh.points()
        arrows = plt_3d.build_arrows(load_matrix, normalized_coord, timing=timing)
        cones = plt_3d.build_cones(restri_matrix, normalized_coord, timing=timing)
    else:
        arrows, cones = None, None
    
    vp_mesh = plt_3d.plot_mesh(mesh, arrows, cones, complete, mesh1)
    
    if len(freq_range) == 3:
        app = pg.mkQApp()
        pg.setConfigOption("background", 'w')
        pg.setConfigOption("foreground", 'k')
        if node_plot is None:
            node_plot = load_matrix[0, :4].reshape(1, 4)
        disp_plot = fc.freqresponse(stif_matrix, mass_matrix, load_vector, free_ind, ngl, alpha, beta, eta, freq_range, node_plot)
        win, freq_plot = plt_3d.plot_freqresponse(freq_range, disp_plot.real)

    # Save
    # if save:
    #     folder_name = "images_3d"
    #     directory = os.path.join(os.path.dirname(__file__), folder_name)
    #     os.makedirs(directory, exist_ok=True)

    #     if len(freq_range) == 3:
    #         save_freq = pg.exporters.ImageExporter(freq_plot)
    #         save_freq.export(os.path.join(directory, "freq.png"))

    #     vp_mesh.screenshot(os.path.join(directory,"mesh_3d.png"))

    # Plot
    vp_mesh.show(viewup="z", axes=4)
    if len(freq_range) == 3: 
        app.exec_()