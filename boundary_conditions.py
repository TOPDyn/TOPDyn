from scipy import spatial
import numpy as np

class BoundConditions:
    def __init__(self, three_dim, nelx, nely, nelz, coord, load_matrix, constr_matrix):

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

    def get_matrices(self, matrix, force):
        """ Gets the force matrix or the contraint matrix.

        Args:
            matrix (:obj:`list`): List passed by the user. 
            coord (:obj:`numpy.array`): Coordinates of the element.
            force (:obj:`bool`): True if encountering force matrix.

        Returns:
            load_matrix or constr_matrix.
        """
        if force:
            index_by_coord, index_by_col, np_matrix = self.turn_dict_into_np(matrix, force=True)
            nodes_coord, nodes_col = self.get_nodes(np_matrix, index_by_coord, index_by_col)
            nodes_matrix = self.get_node_matrix(nodes_coord, nodes_col, np_matrix, index_by_coord, index_by_col, self.ind_nodes_matrix_load, self.ind_np_matrix_load, self.total_cols_load)
        else:
            index_by_coord, index_by_col, np_matrix = self.turn_dict_into_np(matrix, force=False)
            nodes_coord, nodes_col = self.get_nodes(np_matrix, index_by_coord, index_by_col)
            nodes_matrix = self.get_node_matrix(nodes_coord, nodes_col, np_matrix, index_by_coord, index_by_col, self.ind_nodes_matrix_constr, self.ind_np_matrix_constr, self.total_cols_constr)
        return nodes_matrix.astype(int)

    def turn_dict_into_np(self, dicti_matrix, force):
        """ Transforms dictionaries into numpy array.

        Args:
            dicti_matrix (:obj:`list`): List of dictionaries passed by the user.
            force (:obj:`bool`): True if encountering force matrix.
            
        Returns:
            numpy array.
        """
        index_by_coord = []
        index_by_col   = []
        matrix = []
        for i, dict_row in enumerate(dicti_matrix):
            if not("eps" in dict_row):
                if force:
                    if self.three_dim:
                        aux = [dict_row["x_coord"], dict_row["y_coord"], dict_row["z_coord"], dict_row["x_direc"], dict_row["y_direc"], dict_row["z_direc"], dict_row["force"]]
                    else:
                        aux = [dict_row["x_coord"], dict_row["y_coord"], dict_row["x_direc"], dict_row["y_direc"], dict_row["force"], 0]  
                else:
                    if self.three_dim:
                        aux = [dict_row["x_coord"], dict_row["y_coord"], dict_row["z_coord"], dict_row["constrain_disp_x"], dict_row["constrain_disp_y"], dict_row["constrain_disp_z"]]
                    else:
                        aux = [dict_row["x_coord"], dict_row["y_coord"], dict_row["constrain_disp_x"], dict_row["constrain_disp_y"], 0] 
                index_by_coord.append(i)
            
            elif "eps" in dict_row:
                if force:
                    if self.three_dim:
                        aux = [dict_row["coord"], dict_row["axis"], dict_row["eps"], dict_row["x_direc"], dict_row["y_direc"], dict_row["z_direc"], dict_row["force"]]
                    else:
                        aux = [dict_row["coord"], dict_row["axis"], dict_row["x_direc"], dict_row["y_direc"], dict_row["force"], dict_row["eps"]] 
                else:
                    if self.three_dim:
                        aux = [dict_row["coord"], dict_row["axis"], dict_row["eps"], dict_row["constrain_disp_x"], dict_row["constrain_disp_y"], dict_row["constrain_disp_z"]]
                    else:
                        aux = [dict_row["coord"], dict_row["axis"], dict_row["constrain_disp_x"], dict_row["constrain_disp_y"], dict_row["eps"]]
                index_by_col.append(i)
            
            matrix.append(aux)

        matrix_L = np.array(matrix)
        return index_by_coord, index_by_col, matrix_L

    def get_nodes(self, np_matrix, index_by_coord, index_by_col):
        """ Gets nodes by a coordinate or a column.

        Args:
            coord (:obj:`numpy.array`): Coordinates of the element.
            np_matrix (:obj:`numpy.array`): List passed to an array.
            index_by_coord (:obj:`list`): indices of elements passed by coordinates.
            index_by_col (:obj:`list`): indices of elements passed by columns.

        Returns:
            Nodes.
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
            coord (:obj:`numpy.array`): mesh coordinates.
            coord_user (:obj:`numpy.array`): user coordinates.
            
        Returns:
            Nodes of the coordinates provided.
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
        """ Get node numbers that are equal to coord.

        Args:
            coord (:obj:`numpy.array`): mesh coordinates.
            coord_user (:obj:`numpy.array`): coordinates in one direction (x or y).
            eps (:obj:`float`): Acceptable margin of difference.
            column (:obj:`int`): Direction to compare (x or y).
                x = 1 and y = 2.

        Returns:
            Nodes.
        """
        dif = np.abs(self.coord[:, column] - coord_user)
        mask = dif < eps
        return (self.coord[mask, 0]).astype('int')

    def get_node_matrix(self, nodes_coord, nodes_col, np_matrix, index_by_coord, index_by_col, ind1, ind2, total_cols):
        """ Creates the node matrix. TODO: SE CHAMAVA GET_MATRIX E ESTA A MESMA COISA

        Args:
            nodes_coord (:obj:`list`): Nodes passed by coordinate.
            nodes_col (:obj:`list`): Nodes passed by column.
            coord (:obj:`numpy.array`): Coordinates of the element.
            matrix (:obj:`numpy.array`): List passed to an array.
            index_by_coord (:obj:`list`): indices of elements passed by coordinates.
            index_by_col (:obj:`list`): indices of elements passed by columns.
            ind1 (:obj:`int`): Indices of matrix with nodes.
            ind2 (:obj:`int`): Indices of matrix passed by user.
            total_cols (:obj:`int`): Number of columns desired for the matrix.

        Returns:
            matrix with nodes.
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
        """ Delete specific DOFs from all DOFs. TODO MUDAR ESSE DOF PARA NODES
        
        Args: 
            nelx (:obj:`int`): Number of elements on the x-axis.
            nely (:obj:`int`): Number of elements on the y-axis.
            del_dofs (:obj:`numpy.array`): Array with DOFs to be removed.
        
        Returns:
            Array without the DOFs removed.
        """
        if self.three_dim:
            nodes = np.arange((nelx+1) * (nely + 1) * nelz * 3)
        else:
            nodes = np.arange((nelx+1) * (nely + 1) * 2)
        return np.delete(nodes, del_dofs)

    def get_load_vector(self):
        """ Creates the force vector.

        Args:
            nelx (:obj:`int`): Number of elements on the x-axis.
            nely (:obj:`int`): Number of elements on the y-axis.
            nelz (:obj:`int`): Number of elements on the z-axis.
            load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value.

        Returns:
            Loading vector.
        """
        self.load_vector = np.zeros(self.ngl, dtype=complex)
        if len(self.load_matrix) > 1:
            aux_load_matrix = self.load_matrix[np.argsort(self.load_matrix[:, 0])]
        else:
            aux_load_matrix = self.load_matrix
        force_ind = self.get_dofs(aux_load_matrix)
        self.load_vector[force_ind] = self._duplicate_load(aux_load_matrix)

    def _duplicate_load(self, aux_load_matrix):
        if self.three_dim:
            return self._duplicate_load_3D(aux_load_matrix)
        else:
            return self._duplicate_load_2D(aux_load_matrix)

    def _duplicate_load_2D(self, aux_load_matrix):
        """ Doubled force value for x direction = y direction = 1.

        Args:
            aux_load_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value.
        
        Returns:
            Load.
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
        """ Doubled force value.

        Args:
            aux_load_matrix (:obj:`numpy.array`): Load.
        
        Returns:
            Load values.
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