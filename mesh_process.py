import numpy as np

class Mesh:
    """ Definition of boundary conditions """
    def __init__(self, three_dim, nelx, nely, nelz, lx, ly, lz):
        """
        Args:
            three_dim (:obj:`bool`): True if it's used three dimensional mesh.
            nelx (:obj:`int`): Number of elements on the X-axis. 
            nely (:obj:`int`): Number of elements on the Y-axis. 
            nelz (:obj:`int`): Number of elements on the Z-axis.
            lx (:obj:`float`): X-axis length.
            ly (:obj:`float`): Y-axis length.
            lz (:obj:`float`): Z-axis length.
        """
        self.three_dim = three_dim
        self.nelx = nelx
        self.nely = nely
        self.lx = lx
        self.ly = ly

        if self.three_dim:
            self.nelz = nelz
            self.lz = lz
            self.regular_mesh_H8()
        else:
            self.regular_mesh_Q4()

        if self.three_dim:
            self.dofs = 3
            self.edofs = 24
        else:
            self.dofs = 2
            self.edofs = 8

        self.ind_dofs = self._get_ind_dofs()
        self.generate_ind_rows_cols()

    def generate_ind_rows_cols(self):
        """ Node indexes to make the assembly. """
        vect_indices = self.ind_dofs.flatten()
        self.ind_rows = ((np.tile(vect_indices, (self.edofs,1))).T).flatten()
        self.ind_cols = (np.tile(self.ind_dofs, self.edofs)).flatten()

    def _get_ind_dofs(self):
        """ Get DOFs of each node.
        
        Returns:
            ind_dofs (:obj:`numpy.array`): Nodes DOFs.
        """
        if self.three_dim:
            ind_dofs = (np.array([self.dofs*self.connect[:,1]-1, self.dofs*self.connect[:,1], self.dofs*self.connect[:,1]+1,
                                self.dofs*self.connect[:,2]-1, self.dofs*self.connect[:,2], self.dofs*self.connect[:,2]+1,
                                self.dofs*self.connect[:,3]-1, self.dofs*self.connect[:,3], self.dofs*self.connect[:,3]+1,
                                self.dofs*self.connect[:,4]-1, self.dofs*self.connect[:,4], self.dofs*self.connect[:,4]+1,
                                self.dofs*self.connect[:,5]-1, self.dofs*self.connect[:,5], self.dofs*self.connect[:,5]+1,
                                self.dofs*self.connect[:,6]-1, self.dofs*self.connect[:,6], self.dofs*self.connect[:,6]+1,
                                self.dofs*self.connect[:,7]-1, self.dofs*self.connect[:,7], self.dofs*self.connect[:,7]+1,
                                self.dofs*self.connect[:,8]-1, self.dofs*self.connect[:,8], self.dofs*self.connect[:,8]+1], dtype=int)-2).T
        else:
            ind_dofs = (np.array([self.dofs*self.connect[:,1]-1, self.dofs*self.connect[:,1], self.dofs*self.connect[:,2]-1, self.dofs*self.connect[:,2],
                                self.dofs*self.connect[:,3]-1, self.dofs*self.connect[:,3], self.dofs*self.connect[:,4]-1, self.dofs*self.connect[:,4]], dtype=int)-1).T
        return ind_dofs

    def regular_mesh_Q4(self):
        """ Create a regular Q4 mesh: connectivity and coordinates matrices. """
        # processing of nodal coordinates matrix
        x, y = self.generate_xy_coord()
        #x, y = np.arange(0,lx+dx,dx), np.arange(0,ly+dy,dy)
        nx, ny = len(x), len(y)
        mat_x = (x.reshape(nx, 1)@np.ones((1, ny))).T
        mat_y = y.reshape(ny, 1)@np.ones((1, nx))
        x_t, y_t = mat_x.flatten(), mat_y.flatten()
        ind_coord = np.arange(1, nx*ny+1, 1, dtype=int)
        self.coord = (np.array([ind_coord, x_t, y_t])).T
        # processing of connectivity matrix
        ind_connect = np.arange(1, self.nelx*self.nely+1, 1, dtype=int)
        mat_aux = ind_connect.reshape(self.nely, self.nelx)
        a = np.arange(0, self.nely, 1)
        a = (a.reshape(len(a), 1))@np.ones((1, self.nelx))
        b = (mat_aux + a).flatten()
        self.connect = np.array([ind_connect, b, b+1, b+(self.nelx+2), b+(self.nelx+1)], dtype=int).T

    def regular_mesh_H8(self):
        """ Creates a regular H8 mesh connectivity and coordinates matrices. """
        x, y, z = np.linspace(0, self.lx, num=self.nelx + 1), np.linspace(0, self.ly, num=self.nely + 1), np.linspace(0, self.lz, num=self.nelz + 1)
        nx, ny, nz = len(x), len(y), len(z)
        mat_x = (x.reshape(nx, 1)@np.ones((1, ny*nz))).T
        mat_y = y.reshape(ny, 1)@np.ones((1, nx))
        mat_z = z.reshape(nz, 1)@np.ones((1, nx*ny))
        x_t, y_t, z_t = mat_x.flatten(), np.tile(mat_y.flatten(), nz), mat_z.flatten()
        ind_coord = np.arange(1, (nz)* nx * ny + 1, 1, dtype=int)
        self.coord = (np.array([ind_coord, x_t, y_t, z_t])).T
        # processing of connectivity matrix
        ind_connect = np.arange(1, self.nelz * self.nelx * self.nely + 1, dtype=int)
        mat_aux = ind_connect.reshape(self.nely, self.nelx, self.nelz)
        a = np.arange(0, self.nely * self.nelz, 1)
        for ind in range(self.nely, len(a), self.nely):
            a[ind:] += self.nelx + 1
        c = (a.reshape(len(a),1)@np.ones((1, self.nelx))).reshape(self.nely, self.nelx, self.nelz)
        b = (mat_aux + c).flatten()

        self.connect = np.array([ind_connect, b+(self.nelx+1), b, b+1, b+(self.nelx+2), \
                            b+(self.nelx+1)*(self.nely+1)+(self.nelx+1), b+(self.nelx+1)*(self.nely+1), \
                            b+1+(self.nelx+1)*(self.nely+1), b+(self.nelx+1)*(self.nely+1)+(self.nelx+2)], dtype=int).T    

    def generate_xy_coord(self):
        """ Create the mesh coordinate range. 

        Returns:
            A tuple with X-axis and Y-axis ranges.    
        """
        dx, dy = self.lx/self.nelx, self.ly/self.nely
        return dx * np.arange(self.nelx + 1), dy * np.arange(self.nely + 1)