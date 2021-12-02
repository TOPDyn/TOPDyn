import numpy as np
from scipy.sparse import csc_matrix

class Filters():
    def __init__(self, fac_ratio, x_min_k, lx, nelx, nely, coord, connect, dens_filter):
        self.x_min_k = x_min_k
        self.radius = fac_ratio * lx/nelx
        self.nelx = nelx
        self.nely = nely
        self.coord = coord
        self.connect = connect
        self.dens_filter = dens_filter

        self.get_neighbors_radius()

    def get_neighbors_radius(self):
        """ Check neighboring elements that have the centroid within the predetermined radius. """
        el_number = self.nelx * self.nely
        idx = self.connect[:, 1:] - 1 
        self.centroids = np.empty((el_number, 2))  
        self.centroids[:, 0] = np.sum(self.coord[idx, 1], axis = 1)/4
        self.centroids[:, 1] = np.sum(self.coord[idx, 2], axis = 1)/4
        
        ind_rows = []
        ind_cols = []
        data = []
        cols = 0
        neighbors = []
        for el in range(el_number):
            distance = np.sqrt(np.sum((self.centroids[el] - self.centroids)**2, axis=1))
            mask = distance <= self.radius
            neighbor = mask.nonzero()[0] + 1
            neighbors.extend(neighbor - 1)
            
            hi      = self.radius - distance
            hi_max  = np.maximum(0, hi)
            data.extend(hi_max[mask])
            aux     = len(hi_max[mask])
            rows    = np.repeat(el, aux) #.tolist()
            columns = np.arange(0, aux)
            ind_rows.extend(rows) 
            ind_cols.extend(columns)
            
            if aux > cols:
                cols = aux
        self.H = csc_matrix((data, (ind_rows, ind_cols)), shape=(self.nelx*self.nely, cols)).toarray()
        self.neighbors = csc_matrix((neighbors, (ind_rows, ind_cols)), shape=(self.nelx*self.nely, cols), dtype='int').toarray()

    def apply_filter(self, xval, dfdx, df0dx, df0dx2=None):
        if self.dens_filter:
            dfdx, df0dx, df0dx2 = self._density_filter(dfdx, df0dx, df0dx2)
        else:
            dfdx, df0dx, df0dx2 = self._sensitivity_filter(xval, dfdx, df0dx, df0dx2)
        return dfdx, df0dx, df0dx2

    def _density_filter(self, dfdx, df0dx, df0dx2=None):
        """ Applies the density filter to the derivative of the functions (constrain, objective and multiobjective).

        Args:
            dfdx (:obj:`numpy.array`): Value of the constraint derivative.
            df0dx (:obj:`numpy.array`): Derivative value of the objective function.
            df0dx2 (:obj:`numpy.array`, optional): Derivative value of the second objective function.

        Returns:
            Density filter applied to the derivative values.
        """
        new_deriv_f1, deriv_f, cols = self.set_in_deriv(dfdx, df0dx, df0dx2)
        for el in range(deriv_f.shape[0]):
            idx = self.neighbors[el, :]
            hj  = self.H[idx, :].sum(axis=1)
            for ind in range(cols):
                new_deriv_f1[el, ind] = np.sum((1/hj) * self.H[el, :] * deriv_f[idx, ind])
        return self.set_out_deriv(new_deriv_f1, dfdx, df0dx2)

    def _sensitivity_filter(self, xval, dfdx, df0dx, df0dx2=None):
        """ Applies the sensitivity filter to the derivative of the functions (constrain, objective and multiobjective).

        Args:
            TODO: ADICIONAR XVAL
            dfdx (:obj:`numpy.array`): Value of the constraint derivative.
            df0dx (:obj:`numpy.array`): Derivative value of the objective function.
            df0dx2 (:obj:`numpy.array`, optional): Derivative value of the second objective function.

        Returns:
            Sensitivity filter applied to the derivative values.
        """
        xval2 = xval.copy()
        xval2[xval2 <= self.x_min_k] = self.x_min_k
        
        aux1 = self.H * xval2[self.neighbors.flatten()].reshape(self.H.shape)
        aux3 = 1/(np.sum(self.H, axis=1) * xval2[:, 0])

        new_deriv_f, deriv_f, cols = self.set_in_deriv(dfdx, df0dx, df0dx2)
        for col in range(cols):
            aux2 = aux1 * deriv_f[self.neighbors.flatten(), col].reshape(self.H.shape)
            new_deriv_f[:, col] = aux3 * np.sum(aux2, axis=1)
        return self.set_out_deriv(new_deriv_f, dfdx, df0dx2)

    def set_in_deriv(self, dfdx, df0dx, df0dx2):
        """ Auxiliary function. """
        if df0dx2 is not None:
            cols = 2 + dfdx.shape[0]
        else:
            cols = 1 + dfdx.shape[0]
        
        all_deriv_f = np.empty((df0dx.shape[0], cols))
        for i in range(dfdx.shape[0]):
            all_deriv_f[:, i] = dfdx[i, :]
        
        if df0dx2 is not None:
            all_deriv_f[:, cols-1] = df0dx2[:, 0]
            all_deriv_f[:, cols-2] = df0dx[:, 0]
        else:
            all_deriv_f[:, cols-1] = df0dx[:, 0]
        return np.empty((df0dx.shape[0], cols)), all_deriv_f, cols

    def set_out_deriv(self, all_deriv, dfdx, df0dx2):
        """ Auxiliary function. """
        aux = dfdx.shape[0]
        if df0dx2 is None:
            return all_deriv[:, :aux].T, all_deriv[:, aux].reshape(-1, 1), None
        else:
            return all_deriv[:, :aux].T, all_deriv[:, aux].reshape(-1, 1), all_deriv[:, aux+1].reshape(-1, 1)

    # TODO: Talvez colocar aqueles parametros do beam da otimização aqui dentro, aí se fazer assim criar um novo processo de otimização para 
    # o GCMMA no qual adiciona os outros parametros. Tomar cuidado para deixar os parametros od while de fora dessa classe
    # get_neighbors_radius ok
    # calc_area ok
    # calc_xnew ok
    # set_initxval ok
    # get_passive_el ok
    # set_passive_el ok



