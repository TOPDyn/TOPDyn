from scipy import spatial
import numpy as np

class BoundConditions:
    """ Definition of boundary conditions """
    def __init__(self, three_dim, nelx, nely, nelz, coord, load_matrix, constr_matrix):
        """
        Args:
            three_dim (:obj:`bool`): True if it's used three dimensional mesh. It's not available yet.
            nelx (:obj:`int`): Number of elements on the x-axis. 
            nely (:obj:`int`): Number of elements on the y-axis. 
            nelz (:obj:`int`): Number of elements on the z-axis. It's not available yet.
            coord (:obj:`numpy.array`): mesh coordinates. 
            load_matrix (:obj:`numpy.array`): It's a list of lists. The list can be:
                * [x_coordinate, y_coordinate, force_applied_x, force_applied_y, force_value]
                * [value_coordinate, column_to_compare, force_applied_x, force_applied_y, force_value, error_margin]
            constr_matrix (:obj:`numpy.array`, optional): It's a list of lists. Defaults to None. 
                * [x_coordinate, y_coordinate, constrain_disp_x, constrain_disp_y]
                * [value_coordinate, column_to_compare, constrain_disp_x, constrain_disp_y, error_margin]
        """
        self.three_dim = three_dim
        if self.three_dim:
            self.coord_np_matrix = [0,1,2]
            self.ind_eps = 2

            self.ind_nodes_matrix_constr = [1,2,3]
            self.ind_np_matrix_constr = [3,4,5]
            self.total_cols_constr = 4

            self.ind_nodes_matrix_load = [1,2,3,4]
            self.ind_np_matrix_load = [3,4,5,6]
            self.total_cols_load = 5

            self.ngl = 3 * ((nelx + 1) * (nely + 1) * (nelz + 1))
            self.dofs = 3
        else:
            self.coord_np_matrix = [0,1]
            self.ind_eps = -1

            self.ind_nodes_matrix_constr = [1,2]
            self.ind_np_matrix_constr = [2,3]
            self.total_cols_constr = 3

            self.ind_nodes_matrix_load = [1,2,3]
            self.ind_np_matrix_load = [2,3,4]
            self.total_cols_load = 4

            self.ngl = 2 * ((nelx + 1) * (nely + 1))
            self.dofs = 2

        self.coord = coord

        # Load matrix
        self.load_matrix = self.get_matrices(load_matrix, True)
            
        # Get free indexes
        self.free_ind = None
        if constr_matrix is not None:
            self.constr_matrix = self.get_matrices(constr_matrix, False)
            self.restricted_dofs = self.get_dofs(self.constr_matrix)
            self.free_ind = self.remove_dofs(nelx, nely, nelz, self.restricted_dofs)
        
        # Load vector
        self.get_load_vector()

    def get_matrices(self, matrix, load):
        """ Gets the load matrix or the contraint matrix.

        Args:
            matrix (:obj:`list`): List given by the user. 
            load (:obj:`bool`): True if encountering load matrix.

        Returns:
            nodes_matrix (:obj:`Numpy array`): Boundary conditions. It can be load or node displacement constraint.
        """
        index_by_coord = []
        index_by_col = []
        for i, line in enumerate(matrix):
            if np.isnan(line).any():
                index_by_coord.append(i)
            else:
                index_by_col.append(i)
        if load:
            nodes_coord, nodes_col = self.get_nodes(matrix, index_by_coord, index_by_col)
            nodes_matrix = self.get_node_matrix(nodes_coord, nodes_col, matrix, index_by_coord, index_by_col, self.ind_nodes_matrix_load, self.ind_np_matrix_load, self.total_cols_load)
        else:
            nodes_coord, nodes_col = self.get_nodes(matrix, index_by_coord, index_by_col)
            nodes_matrix = self.get_node_matrix(nodes_coord, nodes_col, matrix, index_by_coord, index_by_col, self.ind_nodes_matrix_constr, self.ind_np_matrix_constr, self.total_cols_constr)
        return nodes_matrix.astype(int)

    def get_nodes(self, np_matrix, index_by_coord, index_by_col):
        """ Gets nodes by a coordinate or a column.

        Args:
            coord (:obj:`numpy.array`): Coordinates of the element.
            np_matrix (:obj:`numpy.array`): Boundary condition matrix.
            index_by_coord (:obj:`list`): Indexes of elements given by coordinates.
            index_by_col (:obj:`list`): Indexes of elements given by columns.

        Returns:
            nodes_coord (:obj:`list`): Nodes of elements given by coordinates.
            nodes_col (:obj:`list`): Nodes of elements given by columns.
        """
        nodes_coord = []
        nodes_col = []

        if len(index_by_coord) > 0:   
            nodes_coord = self.get_nodes_by_coord(np_matrix[np.ix_(index_by_coord, self.coord_np_matrix)])

        if len(index_by_col) > 0:
            for index in index_by_col:
                aux = self.get_nodes_by_col(np_matrix[index, 0], np_matrix[index, self.ind_eps], int(np_matrix[index, 1]))
                nodes_col.append(aux)
        return nodes_coord, nodes_col 

    def get_nodes_by_coord(self, coord_user):
        """ Get node number by coordinate.

        Args:
            coord (:obj:`numpy.array`): Mesh coordinates.
            coord_user (:obj:`numpy.array`): Coordinate given by the user.
            
        Returns:
            nodes (:obj:`numpy.array`):Nodes that match the given coordinate.
        """
        if self.three_dim:
            mesh_coord = self.coord[:, [1,2,3]]
        else:
            mesh_coord = self.coord[:, [1,2]]
        mytree = spatial.cKDTree(mesh_coord)
        _, ind_nodes = mytree.query(coord_user)
        nodes = self.coord[ind_nodes, 0]
        return nodes

    def get_nodes_by_col(self, coord_user, eps, column):
        """ Get node numbers that are equal to one coordinate matrix item.

        Args:
            coord (:obj:`numpy.array`): Mesh coordinates.
            coord_user (:obj:`numpy.array`): Coordinate in one direction (x or y).
            eps (:obj:`float`): Acceptable margin of difference.
            column (:obj:`int`): Direction to compare (x or y).
                x = 1 and y = 2.

        Returns:
            nodes (:obj:`numpy.array`): Array with nodes.
        """
        dif = np.abs(self.coord[:, column] - coord_user)
        mask = dif < eps
        return (self.coord[mask, 0]).astype('int')

    def get_node_matrix(self, nodes_coord, nodes_col, np_matrix, index_by_coord, index_by_col, ind1, ind2, total_cols):
        """ Creates node matrix. TODO

        Args:
            nodes_coord (:obj:`list`): Nodes given by coordinate.
            nodes_col (:obj:`list`): Nodes given by column.
            coord (:obj:`numpy.array`): Coordinates.
            np_matrix (:obj:`numpy.array`): List given to an array.
            index_by_coord (:obj:`list`): indices of elements given by coordinates.
            index_by_col (:obj:`list`): indices of elements given by columns.
            ind1 (:obj:`int`): Indices of matrix with nodes.
            ind2 (:obj:`int`): Indices of matrix given by the user.
            total_cols (:obj:`int`): Number of columns desired for the matrix.

        Returns:
            matrix (:obj:`numpy.array`): Load matrix with nodes.
        """
        if len(nodes_col) > 0:
            len_col = sum([len(listElem) for listElem in nodes_col])
        else:
            len_col = 0
        matrix = np.empty((len(nodes_coord) + len_col, total_cols))

        if len(index_by_coord) > 0:
            matrix[:len(nodes_coord), 0] = nodes_coord
            matrix[:len(nodes_coord), ind1] = np_matrix[np.ix_(index_by_coord, ind2)]
        
        if len(index_by_col) > 0:
            aux = 0
            for i, nodes in enumerate(nodes_col):
                matrix[len(nodes_coord)+aux:len(nodes)+aux+len(nodes_coord), 0] = nodes
                matrix[len(nodes_coord)+aux:len(nodes)+aux+len(nodes_coord), ind1] = np_matrix[index_by_col[i], ind2] 
                aux += len(nodes)
        return matrix

    def get_dofs(self, nodes_dir):
        """ Get DOFs that meet the specified direction.

        Args:
            nodes_dir (:obj:`numpy.array`): matrix with node numbers.
                * x direction, y direction and z dicrfection can be -1, 0 or 1.
        
        Returns: 
            DOFs of each node in array nodes.
        """
        all_dofs = []
        mask = abs(nodes_dir[:, 1]) == 1
        all_dofs.extend(self.dofs * (nodes_dir[mask, 0] - 1))

        mask = abs(nodes_dir[:, 2]) == 1
        all_dofs.extend((self.dofs * nodes_dir[mask, 0]) - 1)

        if self.three_dim:
            mask = abs(nodes_dir[:, 3]) == 1
            all_dofs.extend((self.dofs * nodes_dir[mask, 0]) - 1)

        all_dofs.sort()
        return np.array(all_dofs, dtype='int')
    
    def remove_dofs(self, nelx, nely, nelz, del_dofs):
        """ Delete specific nodes from all nodes.
        
        Args: 
            nelx (:obj:`int`): Number of elements on the x-axis.
            nely (:obj:`int`): Number of elements on the y-axis.
            del_nodes (:obj:`numpy.array`): Nodes to be removed.
        
        Returns:
            Array without the nodes removed.
        """
        if self.three_dim:
            nodes = np.arange((nelx+1) * (nely + 1) * nelz * 3)
        else:
            nodes = np.arange((nelx+1) * (nely + 1) * 2)
        return np.delete(nodes, del_dofs)

    def get_load_vector(self):
        """ Creates load vector. """
        self.load_vector = np.zeros(self.ngl, dtype=complex)
        if len(self.load_matrix) > 1:
            aux_load_matrix = self.load_matrix[np.argsort(self.load_matrix[:, 0])]
        else:
            aux_load_matrix = self.load_matrix
        force_ind = self.get_dofs(aux_load_matrix)
        self.load_vector[force_ind] = self._duplicate_load(aux_load_matrix)

    def _duplicate_load(self, aux_load_matrix):
        """ Duplicate load value.

        Args:
            aux_load_matrix (:obj:`numpy.array`): It can be a 2D or 3D matrix.
        Returns:
            Load vector.
        """
        if self.three_dim:
            return self._duplicate_load_3D(aux_load_matrix)
        else:
            return self._duplicate_load_2D(aux_load_matrix)

    def _duplicate_load_2D(self, aux_load_matrix):
        """ Duplicate load value for x direction = y direction = 1.

        Args:
            aux_load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value.
        
        Returns:
            load_values (:obj:`numpy.array`): Load vector.
        """
        mask = ((abs(aux_load_matrix[:, 1]) == 1.) & (abs(aux_load_matrix[:, 2]) == 1.)).nonzero()[0]
        force_values = aux_load_matrix[:, 3]
        if mask.size != 0:
            force_values = np.insert(force_values, mask, aux_load_matrix[mask, 3])
        # Change force sign
        aux = aux_load_matrix[:, [1,2]].ravel()
        aux = aux[aux!=0]
        if (aux<0).any():
            force_values[aux<0] *= -1
        return force_values

    def _duplicate_load_3D(self, aux_load_matrix):
        """ Doubled load value for x direction = y direction = 1 or
        x direction = y direction = z direction = 1.

        Args:
            aux_load_matrix (:obj:`numpy.array`): Load.
        
        Returns:
            load_vector (:obj:`numpy.array`): Load vector.
        """
        load_vector = aux_load_matrix[:, -1]
        aux = np.sum(abs(aux_load_matrix[:, 1:-1]), axis=1)
        idx_dup = np.where((aux>1) & (aux<3))[0]
        idx_trip = np.repeat(np.where(aux>2), 2)
        if idx_dup.size != 0 or idx_trip.size != 0:
            load_vector = np.insert(load_vector, np.hstack([idx_dup, idx_trip]), np.hstack([aux_load_matrix[idx_dup, -1], aux_load_matrix[idx_trip, -1]]))
        # Change force sign
        aux = aux_load_matrix[:, [1, 2, 3]].ravel()
        aux = aux[aux!=0]
        if (aux<0).any():
            load_vector[aux < 0] *= -1
        return load_vector