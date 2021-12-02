import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as cl

class PlotsFem2d():
    def __init__(self, coord, connect):

        self.coord = coord
        self.connect = connect
       
    def show_nodes(self):
        """ Plot nodes of the mesh. TODO:ESSA FUNÇÃO SERA ALGUMA VEZ USADA? """
        ax = plt.axes()
        ax.scatter(self.coord[:,1], self.coord[:,2])

        for i in range(self.coord.shape[0]):
            ax.annotate(self.coord[i,0], (self.coord[i,1], self.coord[i,2]))

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.show()

    def change_disp_shape(self, disp_vector):
        """ Transform displacement vector in matrix.
        
        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
        
        Returns: 
            Displacement of each axis.
        """
        new_U = np.empty((int(len(disp_vector)/2), 2))
        new_U[:,0] = disp_vector[::2]
        new_U[:,1] = disp_vector[1::2] 
        return new_U

    def apply_disp(self, disp_vector, factor):
        """ Apply displacement to coordinates. 
        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
            coord (:obj:`numpy.array`): mesh coordinates.
            factor (:obj:`float`): Factor to deform the mesh.
        
        Returns: 
            Displaced mesh coordinates.
        """
        new_coord = self.coord.copy()
        new_coord[:, 1:] += disp_vector * factor
        return new_coord

    def build_collection(self, coord_def, disp_vector=None):
        """ Build quad mesh.
        
        Args:
            TODO
            disp_vector (:obj:`numpy.array`, optional): Displacement.
        
        Returns:
            Matplotlib collection object.
        """
        x, y = coord_def[:, 1], coord_def[:, 2]
        xy = np.c_[x, y]
        squares = np.asarray(self.connect[:, 1:])
        verts = xy[squares - 1]
        pc = cl.PolyCollection(verts)
        if disp_vector is not None:
            pc.set_array(disp_vector)
            pc.set_cmap("viridis")
        else:
            pc.set_edgecolor("black")
            pc.set_facecolor("None")
        return pc

    def plot_collection(self, ax, lx, ly, coord_def, pc, load_matrix=None, restri_matrix=None):
        """ Plot mesh, force arrows and constrain nodes. 
        
        Args:
            lx (:obj:`float`): X-axis length.
            ly (:obj:`float`): Y-axis length.
            coord (:obj:`numpy.array`): Coordinates of the element.
            pc (matplotlib.collections): Matplotlib collection object.
            load_matrix (:obj:`numpy.array`, optional): The columns are respectively node, x direction, y direction, force value. 
            restri_matrix (:obj:`numpy.array`, optional)= The columns are respectively node, x direction, y direction. Defaults to None. 
        
        Returns:
            A figure object and a single Axes object from matplotlib.pyplot.
        """
        ax.add_collection(pc)
        ax.autoscale()
        ax.set_aspect('equal')
        x, y = coord_def[:, 1], coord_def[:, 2]
        ax.plot(x, y, ls = "", color = "black")
        max_size, min_size = self._get_size(lx, ly)
        size = 0.06 * coord_def[:, max_size].max() - coord_def[:, max_size].min()
        
        # Plot load arrows
        if load_matrix is not None:
            ax = self._plot_load(ax, coord_def, load_matrix, min_size, size)
        
        # Plot restrict nodes 
        if restri_matrix is not None:
            ax = self._plot_constr_nodes(ax, coord_def, restri_matrix)
        
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_xlim(coord_def[:, 1].min() - size, coord_def[:, 1].max() + size)
        ax.set_ylim(coord_def[:, 2].min() - size, coord_def[:, 2].max() + size)   

    def _get_size(self, lx, ly):
        """ Gets columns with maximum and minimum length.

        Args:
            lx (:obj:`float`): X-axis length.
            ly (:obj:`float`): Y-axis length.
        Returns:
            Column indexes with maximum and minimum length, respectively.
        """
        if lx > ly:
            max_size = 1
            min_size = 2
        else:
            max_size = 2
            min_size = 1
        return max_size, min_size

    def _plot_load(self, ax, coord, load_matrix, min_size, size):
        """ Add load vector to the plot of deformed mesh.
            
        Args:    
            ax (matplotlib.axes.Axes): Matplotlib axes object.
            coord (:obj:`numpy.array`): mesh coordinates.
            load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value. 
            min_size (:obj:`float`): Value to define the width of load arrow.
            size (:obj:`float`): Load arrow length.
        
        Returns:
            Matplotlib axes object.
        """
        factor = coord[:, min_size].max() - coord[:, min_size].min()
        ind = (load_matrix[:, 0] - 1).astype('int')
        for i in range(load_matrix.shape[0]): 
            if load_matrix[i, 1] == 1:
                ax.arrow(coord[ind[i], 1], coord[ind[i], 2], size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
            if load_matrix[i, 1] == -1:
                ax.arrow(coord[ind[i], 1], coord[ind[i], 2], -size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
            
            if load_matrix[i, 2] == -1:
                ax.arrow(coord[ind[i], 1], coord[ind[i], 2], 0, -size, shape='full', length_includes_head=True, width = factor * 0.01)
            
            if load_matrix[i, 2] == 1:
                ax.arrow(coord[ind[i], 1], coord[ind[i], 2], 0, size, shape='full', length_includes_head=True, width = factor * 0.01)
        return ax

    def _plot_constr_nodes(self, ax, coord, restri_matrix):
        """ Add contrain nodes to the plot of deformed mesh.
            
        Args:    
            ax (matplotlib.axes.Axes): Matplotlib axes object.
            coord (:obj:`numpy.array`): mesh coordinates.
            restri_matrix (numpy.array)= The columns are respectively node, x direction, y direction.
        
        Returns:
            Matplotlib axes object.
        """
        both_restri = (restri_matrix[:, 1] == 1) & (restri_matrix[:, 2] == 1)
        ind = restri_matrix[both_restri, 0] - 1
        ax.scatter(coord[ind, 1], coord[ind, 2], marker=(3, 0, 0), s=120, color = 'red', linestyle='None')
        dist = 0.01 * (coord[1,0] - coord[0,0])
        ax.scatter(coord[ind, 1] - dist, coord[ind, 2], marker=(3, 0, 270), s=120, color = 'green', linestyle='None')
        
        x_restri = (restri_matrix[:, 1] == 1) & (restri_matrix[:, 2] == 0)
        ind = restri_matrix[x_restri, 0] - 1
        ax.scatter(coord[ind, 1], coord[ind, 2], marker=(3, 0, 270), s=120, color = 'green', linestyle='None')
        
        y_restri = (restri_matrix[:, 1] == 0) & (restri_matrix[:, 2] == 1)
        ind = restri_matrix[y_restri, 0] - 1
        plt.scatter(coord[ind, 1], coord[ind, 2], marker=(3, 0, 0), s=120, color = 'red', linestyle='None')
        return ax

    def plot_freq_rsp(self, ax, node, freq_range, disp_vector):
        """ Plot the frequency response.
                
        Args:    
            node (:obj:`list`): Node that was calculated the frequency response.
            freq_range (:obj:`list`): Range of analyzed frequencies.
                
                * First value is the minimum frequency.
                * Second value is the maximum frequency.
                * Third value is the step between each calculation of the function. 
            disp_vector (:obj:`numpy.array`): Displacement.
        
        Returns:
            A figure object and a single Axes object from matplotlib.pyplot.
        """
        x = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
        y = 10 * np.log10(abs(disp_vector))
        ax.plot(x, y)
        ax.set_title('Node ' + str(node[0,0]), fontsize=15)
        ax.set_xlabel('frequency [Hz]', fontsize=16)
        ax.set_ylabel('displacement [N]', fontsize=16)
        ax.set_xlim(0, x[-1])
        ax.grid()

    def save_figs(self, mesh, freq_rsp):
        """ Saves figure.

        Args:
            fig: Object with the graph. 
        """   
        folder_name = 'images'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        os.makedirs(directory, exist_ok=True)
        mesh.savefig(os.path.join(directory, "mesh.png"))
        freq_rsp.savefig(os.path.join(directory, "freq_rsp.png"))