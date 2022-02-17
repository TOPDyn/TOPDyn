import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as cl

class PlotsFem2d():
    def __init__(self, coord, connect):

        self.coord = coord
        self.connect = connect

       
    def show_nodes(self):
        """ Plot nodes of the mesh. """
        ax = plt.axes()
        ax.scatter(self.coord[:,1], self.coord[:,2])

        for i in range(self.coord.shape[0]):
            ax.annotate(self.coord[i,0], (self.coord[i,1], self.coord[i,2]))

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.show()

    def change_disp_shape(self, disp_vector):
        """ Transform 1D displacement vector in 2D.
        
        Args:
            disp_vector (:obj:`numpy.array`): 1D displacement.
        
        Returns: 
            new_U (:obj:`numpy.array`) 2D displacement.
        """
        new_U = np.empty((int(len(disp_vector)/2), 2))
        new_U[:,0] = disp_vector[::2]
        new_U[:,1] = disp_vector[1::2] 
        return new_U

    def apply_disp(self, disp_vector, factor):
        """ Apply displacement to coordinates. 

        Args:
            disp_vector (:obj:`numpy.array`): Displacement.
            factor (:obj:`float`): Factor to deform the mesh.
        """
        self.coord_U = self.coord.copy()
        self.coord_U[:, 1:] += disp_vector * factor

    def remove_nodes(self, passive_coordinate):
        """ Get passive elements and create matrix without passive elements. 

        Args:
            passive_coordinates (:obj:`tuple`, optional): Region that the shape will not be changed.
        """
        condition1 = np.all([passive_coordinate[0][0] <= self.coord[:,1], self.coord[:,1] <= passive_coordinate[0][1]], axis=0)
        condition2 = np.all([passive_coordinate[1][0] <= self.coord[:,2], self.coord[:,2] <= passive_coordinate[1][1]], axis=0)
        condition = np.all([condition1, condition2], axis=0)

        # WHY NOT?
        # condition = (passive_coordinate[0][0] <= coord[:,1] and coord[:,1] <= passive_coordinate[0][1]) or \
        #             (passive_coordinate[1][0] <= coord[:,2] and coord[:,2] <= passive_coordinate[1][0])

        self.passive_nodes = self.coord[condition, 0]
        ind_cond = []
        for i, row in enumerate(self.connect[:, 1:]):
            if not (np.isin(row, self.passive_nodes).any()): #np.in1d(row, nodes).any() FASTER?
                ind_cond.append(i)
        self.passive_connect = self.connect[ind_cond, :]

    def build_collection(self, passive_coordinates=None, disp_vector=None):
        """ Build quad mesh.
        
        Args:
            passive_coordinates (:obj:`tuple`, optional): Region that the shape will not be changed.
            disp_vector (:obj:`numpy.array`, optional): 1D Displacement.
        """
        if passive_coordinates is not None:
            self.remove_nodes(passive_coordinates)
            verts = np.empty((len(self.passive_connect), 4, 2))
            for i, nodes in enumerate(self.passive_connect):
                node_ordered = np.empty(4, dtype=int)
                node_ordered[0] = np.where(self.coord_U[:, 0] == nodes[1])[0]
                node_ordered[1] = np.where(self.coord_U[:, 0] == nodes[2])[0]
                node_ordered[2] = np.where(self.coord_U[:, 0] == nodes[3])[0]
                node_ordered[3] = np.where(self.coord_U[:, 0] == nodes[4])[0]
                verts[i, :, :] = self.coord_U[node_ordered, 1:]  
        else:
            verts = self.coord_U[:, 1:][self.connect[:, 1:] - 1]
        pc = cl.PolyCollection(verts)
        if disp_vector is not None:
            pc.set_array(disp_vector)
            pc.set_cmap("viridis")
        else:
            pc.set_edgecolor("black")
            pc.set_facecolor("None")
        self.collection = pc

    def plot_collection(self, ax, lx, ly, load_matrix=None, restri_matrix=None, passive_el=False):
        """ Plot mesh, force arrows and constrain nodes. 
        
        Args:
            ax (:obj:`matplotlib.axes.Axes`): Axis on which the collection will be plotted.
            lx (:obj:`float`): X-axis length.
            ly (:obj:`float`): Y-axis length.
            load_matrix (:obj:`numpy.array`, optional): The columns are respectively node, x direction, y direction, force value. 
            restri_matrix (:obj:`numpy.array`, optional): The columns are respectively node, x direction, y direction. 
            passive_el (:obj:`boolean`, optional): If True the mesh contain passive elements.
        """
        ax.add_collection(self.collection)
        ax.autoscale()
        ax.set_aspect('equal')
        x, y = self.coord_U[:, 1], self.coord_U[:, 2]
        ax.plot(x, y, ls = "", color = "black")
        max_size, min_size = self._get_size(lx, ly)
        size = 0.06 * self.coord_U[:, max_size].max() - self.coord_U[:, max_size].min()
        nodes = self.passive_nodes if passive_el else None

        # Plot load arrows
        if load_matrix is not None:
            ax = self._plot_load(ax, load_matrix, min_size, size, nodes)
        
        # Plot restrict nodes 
        if restri_matrix is not None:
            ax = self._plot_constr_nodes(ax, restri_matrix, nodes)
        
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_xlim(self.coord_U[:, 1].min() - size, self.coord_U[:, 1].max() + size)
        ax.set_ylim(self.coord_U[:, 2].min() - size, self.coord_U[:, 2].max() + size)   

    def _get_size(self, lx, ly):
        """ Gets columns with maximum and minimum length.

        Args:
            lx (:obj:`float`): X-axis length.
            ly (:obj:`float`): Y-axis length.

        Returns:
            max_size (:obj:`float`): Column index with maximum length.
            min_size (:obj:`float`): Column index with minimum length.
        """
        if lx > ly:
            max_size = 1
            min_size = 2
        else:
            max_size = 2
            min_size = 1
        return max_size, min_size

    def _plot_load(self, ax, load_matrix, min_size, size, passive_nodes):
        """ Add load vector to the plot of deformed mesh.
            
        Args:    
            ax (matplotlib.axes.Axes): Axis on which the collection will be plotted.
            load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value. 
            min_size (:obj:`float`): Value to define the width of load arrow.
            size (:obj:`float`): Load arrow length.
        
        Returns:
            ax (matplotlib.axes.Axes): Axis on which the collection will be plotted.
        """
        factor = self.coord_U[:, min_size].max() - self.coord_U[:, min_size].min()
        for load in load_matrix:
            if passive_nodes is not None:
                aux = load[0] not in passive_nodes
                ind = int(load[0] - 1) if aux else False
            else:
                ind = int(load[0] - 1)
            
            if ind:
                if load[1] == 1:
                    ax.arrow(self.coord_U[ind, 1], self.coord_U[ind, 2], size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
                if load[1] == -1:
                    ax.arrow(self.coord_U[ind, 1], self.coord_U[ind, 2], -size, 0, shape='full', length_includes_head=True, width = factor * 0.01)
                
                if load[2] == -1:
                    ax.arrow(self.coord_U[ind, 1], self.coord_U[ind, 2], 0, -size, shape='full', length_includes_head=True, width = factor * 0.01)
                
                if load[2] == 1:
                    ax.arrow(self.coord_U[ind, 1], self.coord_U[ind, 2], 0, size, shape='full', length_includes_head=True, width = factor * 0.01)
        return ax

    def get_constr_nodes(self, node_constr_matrix, type_constrain, passive_nodes):
        """ Get nodes with displacement constraint that doesn't have passive elements. 
        
        Args:
            node_constr_matrix (:obj:`numpy.array`): Constrained displacement nodes.
            type_constrain (:obj:`numpy.array`): Booleand array.
            passive_nodes (:obj:`numpy.array`): Nodes that are in passive region.
        
        Returns:
            ind (:obj:`numpy.array`): Non-passive node indices.
        """
        if passive_nodes is not None:
            matrix_nodes = node_constr_matrix[type_constrain, 0]
            ind = []
            for node in matrix_nodes:
                aux = node not in passive_nodes
                if aux:
                    ind.append(node - 1)
        else:
            ind = node_constr_matrix[type_constrain, 0] - 1
        return ind

    def _plot_constr_nodes(self, ax, node_constr_matrix, passive_nodes):
        """ Add contrain nodes to the plot of deformed mesh.
            
        Args:    
            ax (matplotlib.axes.Axes): Axis on which the collection will be plotted.
            node_constr_matrix (:obj:`numpy.array`): Constrained displacement nodes.
            passive_nodes (:obj:`numpy.array`): Nodes that are in passive region.
        
        Returns:
            ax (matplotlib.axes.Axes): Axis on which the collection will be plotted.
        """
        both_constrained = (node_constr_matrix[:, 1] == 1) & (node_constr_matrix[:, 2] == 1)
        if both_constrained.any():
            ind = self.get_constr_nodes(node_constr_matrix, both_constrained, passive_nodes)
            if len(ind) > 0:
                ax.scatter(self.coord_U[ind, 1], self.coord_U[ind, 2], marker=(3, 0, 0), s=120, color = 'red', linestyle='None')
                dist = 0.01 * (self.coord_U[1,0] - self.coord_U[0,0])
                ax.scatter(self.coord_U[ind, 1] - dist, self.coord_U[ind, 2], marker=(3, 0, 270), s=120, color = 'green', linestyle='None')
        
        x_constrained = (node_constr_matrix[:, 1] == 1) & (node_constr_matrix[:, 2] == 0)
        if x_constrained.any():
            ind = ind = self.get_constr_nodes(self, node_constr_matrix, x_constrained, passive_nodes)
            if len(ind) > 0:
                ax.scatter(self.coord_U[ind, 1], self.coord_U[ind, 2], marker=(3, 0, 270), s=120, color = 'green', linestyle='None')
        
        y_constrained = (node_constr_matrix[:, 1] == 0) & (node_constr_matrix[:, 2] == 1)
        if y_constrained.any():
            ind = ind = self.get_constr_nodes(self, node_constr_matrix, x_constrained, passive_nodes)
            if len(ind) > 0:
                ax.scatter(self.coord_U[ind, 1], self.coord_U[ind, 2], marker=(3, 0, 0), s=120, color = 'red', linestyle='None')
        return ax

    def plot_freq_rsp(self, ax, node, freq_range, disp_vector):
        """ Plot the frequency response.
                
        Args:
            ax (matplotlib.axes.Axes): Axis on which the collection will be plotted.
            node (:obj:`list`): Node that was calculated the frequency response.
            freq_range (:obj:`list`): Range of analyzed frequencies.
                
                * First value is the minimum frequency.
                * Second value is the maximum frequency.
                * Third value is the step between each calculation of the function. 
            disp_vector (:obj:`numpy.array`): Displacement.
        
        Returns:
            ax (matplotlib.axes.Axes): Axis on which the collection will be plotted.
        """
        x = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
        y = 10 * np.log10(abs(disp_vector))
        ax.plot(x, y)
        ax.set_title('Node ' + str(node[0,0]), fontsize=15)
        ax.set_xlabel('frequency [Hz]', fontsize=16)
        ax.set_ylabel('displacement [N]', fontsize=16)
        ax.set_xlim(0, x[-1])
        ax.grid()
        return ax

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