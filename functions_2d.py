from time import time
import numpy as np
import cmath
from scipy import spatial
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

def shapeQ4(ssx, ttx):
    """ Linear Shape Functions and Derivatives.

    Args:
        ssx (:obj:`float`): Local coordinate of the element on the X-axis.
        ttx (:obj:`float`): Local coordinate of the element on the Y-axis.
    
    Returns:
        Tuple with the shape function and its derivative.
    """
    numerator = 4
    #shape functions
    phi = np.array([(1-ssx)*(1-ttx), (1+ssx)*(1-ttx), (1+ssx)*(1+ttx), (1-ssx)*(1+ttx)], dtype=float)/numerator
    #derivatives
    dphi = np.array([[-(1-ttx),  (1-ttx), (1+ttx), -(1+ttx)],
                     [-(1-ssx), -(1+ssx), (1+ssx),  (1-ssx)]], dtype=float)/numerator
    #
    return phi, dphi

def matricesQ4(ee, coord, connect, E, v, rho):
    """ Q4 stiffness and mass matrices.

    Args:
        ee (:obj:`int`): Element.
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.
        rho (:obj:`float`): Density.

    Returns:
        Tuple with stiffness and mass matrices.
    """
    # constitutive info
    tempc = E/(1 - v**2)
    tempn = tempc*(1 - v)/2
    tempt = tempc*v
    cttv = np.array([[tempc, tempt, 0], [tempt, tempc, 0], [0, 0, tempn]], dtype=float)
    # integration points
    nint, con, wps = 4, 1/np.sqrt(3), 1
    pint = np.array([[-con, -con], [con, -con], [con, con], [-con, con]], dtype=float)
    # preallocating elementary matrices
    Ke, Me = 0, 0
    # integration
    for i in range(nint):
        ssx, ttx = pint[i, 0], pint[i, 1]
        phi, dphi = shapeQ4(ssx,ttx)
        ie = connect[ee,1:]-1
        dxdy = dphi@coord[ie, 1:3] # note: dxds, dyds, dxdt, dydt = dxdy[0,0], dxdy[0,1], dxdy[1,0],dxdy[1,1] 
        JAC = np.array([[dxdy[0,0], dxdy[0,1]],[dxdy[1,0], dxdy[1,1]]], dtype=float) # JAC = np.array([[dxds, dyds],[dxdt, dydt]], dtype=float)
        detJAC = JAC[0,0]*JAC[1,1] - JAC[0,1]*JAC[1,0]
        iJAC = (1/detJAC)*np.array([[JAC[1,1], -JAC[0,1]], [-JAC[1,0], JAC[0,0]]], dtype=float)
        dphi_t = iJAC @ dphi
        #
        B = np.array([[dphi_t[0,0], 0, dphi_t[0,1], 0, dphi_t[0,2], 0, dphi_t[0,3], 0],
                      [0, dphi_t[1,0], 0, dphi_t[1,1], 0, dphi_t[1,2], 0, dphi_t[1,3]],
                      [dphi_t[1,0], dphi_t[0,0], dphi_t[1,1], dphi_t[0,1],dphi_t[1,2], dphi_t[0,2], dphi_t[1,3], dphi_t[0,3]]], dtype=float)
        #
        N = np.array([[phi[0], 0, phi[1], 0, phi[2], 0, phi[3], 0],
                      [0, phi[0], 0, phi[1], 0, phi[2], 0, phi[1]]], dtype=float)
        #
        Ke += B.T@(cttv@B)*(detJAC*wps)
        Me += rho*N.T@N*(detJAC*wps)
        #
    return Ke.real, Me.real

def generate_xy_coord(lx, ly, nelx, nely):
    dx, dy = lx/nelx, ly/nely
    return dx * np.arange(nelx + 1), dy * np.arange(nely + 1)

def regularmeshQ4(lx, ly, nelx, nely, timing=False):
    """ Create a regular Q4 mesh.

    Args:
        lx (:obj:`float`): X-axis length.
        ly (:obj:`float`): Y-axis length.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.
    
    Returns:
        Tuple with the coordinate matrix, connectivity, and the indexes of each node.    
    """
    # processing of nodal coordinates matrix
    t0 = time()
    dofs, edofs = 2, 8
    x, y = generate_xy_coord(lx, ly, nelx, nely)
    #x, y = np.arange(0,lx+dx,dx), np.arange(0,ly+dy,dy)
    nx, ny = len(x), len(y)
    mat_x = (x.reshape(nx, 1)@np.ones((1, ny))).T
    mat_y = y.reshape(ny, 1)@np.ones((1, nx))
    x_t, y_t = mat_x.flatten(), mat_y.flatten()
    ind_coord = np.arange(1, nx*ny+1, 1, dtype=int)
    coord = (np.array([ind_coord, x_t, y_t])).T
    # processing of connectivity matrix
    ind_connect = np.arange(1 ,nelx*nely+1, 1, dtype=int)
    mat_aux = ind_connect.reshape(nely, nelx)
    a = np.arange(0, nely, 1)
    a = (a.reshape(len(a), 1))@np.ones((1, nelx))
    b = (mat_aux + a).flatten()
    connect = np.array([ind_connect, b, b+1, b+(nelx+2), b+(nelx+1)], dtype=int).T
    # processing the dofs indices (rows and columns) for assembly
    ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,2]-1, dofs*connect[:,2],
                          dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,4]-1, dofs*connect[:,4]], dtype=int)-1).T
    vect_indices = ind_dofs.flatten()
    ind_rows = ((np.tile(vect_indices, (edofs,1))).T).flatten()
    ind_cols = (np.tile(ind_dofs, edofs)).flatten()
    tf = time()
    if timing:
        print("Time to process mesh: " + str(round((tf - t0),6)) + '[s]')
    #
    return coord, connect, ind_rows, ind_cols
  
def solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, freq, timing=False, **kwargs):
    """ Assembly and solution.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        ind_rows (:obj:`numpy.array`): Node indexes to make the assembly.
        ind_cols (:obj:`numpy.array`): Node indexes to make the assembly.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.  
        rho (:obj:`float`): Density.  
        alpha (:obj:`float`): Damping coefficient proportional to mass. 
        beta (:obj:`float`): Damping coefficient proportional to stiffness.  
        eta (:obj:`float`): Damping coefficient. 
        freq (:obj:`int`): Analyzed frequency.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.

    Returns:
        Displacement array.
    """
    #
    t01 = time()
    ngl = 2 * ((nelx + 1) * (nely + 1))
    data_k = np.zeros((nelx * nely, 64), dtype=float)
    data_m = np.zeros((nelx * nely, 64), dtype=float)
    #
    for el in range(nelx * nely):
        Ke, Me = matricesQ4(el, coord, connect, E, v, rho)
        data_k[el,:] = Ke.flatten() 
        data_m[el,:] = Me.flatten() 
    #
    data_k = data_k.flatten()
    data_m = data_m.flatten()
    stif_matrix = csc_matrix((data_k, (ind_rows, ind_cols)), shape=(ngl, ngl))
    mass_matrix = csc_matrix((data_m, (ind_rows, ind_cols)), shape=(ngl, ngl))
    tf1 = time()
    if timing:
        print("Time to assembly global matrices: " + str(round((tf1 - t01), 6)) + '[s]')
    # harmonic solution: direct method and no damping
    t02 = time()
    omega = 2 * np.pi * freq
    damp_matrix = alpha * mass_matrix + (beta + eta/omega)*stif_matrix
    if kwargs.get('load_vector') is not None:
        load_vector = kwargs.get('load_vector')
    #
    if kwargs.get('unrestricted_ind') is not None:
        free_ind = kwargs.get('unrestricted_ind')
        F = load_vector[free_ind]
        K = stif_matrix[free_ind, :][:, free_ind]
        M = mass_matrix[free_ind, :][:, free_ind]
        C = damp_matrix[free_ind, :][:, free_ind]
        #Kd = -(omega**2)*M + K
        Kd = -(omega**2)*M + 1j*omega*C + K
        U = spsolve(Kd,F)
        disp_vector = np.zeros((ngl), dtype=complex)
        disp_vector[free_ind] = U
    else:
        Kd = -(omega**2)*mass_matrix + 1j*omega*damp_matrix + stif_matrix
        disp_vector = spsolve(Kd,load_vector)

    tf2 = time()
    if timing:
        print("Time to solve the harmonic analysis problem: " + str(round((tf2 - t02), 6)) + '[s]')
        print("Total time elapsed in solution: " + str(round((tf2 - t01), 6)) + '[s]')
    #
    return disp_vector   

def get_nodes_by_coord(coord, coord_user):
    ''' Get node number by coordinate.

    Args:
        coord (:obj:`numpy.array`): mesh coordinates.
        coord_user (:obj:`numpy.array`): user coordinates.
        
    Returns:
        Nodes of the coordinates provided.
    '''
    mytree = spatial.cKDTree(coord[:, [1,2]])
    _, ind_nodes = mytree.query(coord_user)
    nodes = coord[ind_nodes, 0]
    return nodes

def get_nodes1d(coord, coord_user, eps, column):
    ''' Get node numbers that are equal to coord.

    Args:
        coord (:obj:`numpy.array`): mesh coordinates.
        coord_user (:obj:`numpy.array`): coordinates in one direction (x or y).
        eps (:obj:`float`): Acceptable margin of difference.
        column (:obj:`int`): Direction to compare (x or y).
            x = 1 and y = 2.

    Returns:
        Nodes.
    '''
    dif = np.abs(coord[:, column] - coord_user)
    mask = dif < eps

    return (coord[mask, 0]).astype('int')

def get_dofs(nodes_dir):
    ''' Get DOFs that meet the specified direction.

    Args:
        nodes_dir (:obj:`numpy.array`): [nodes numbers, x_direction, y_direction].
            x_direction and y_direction can be -1, 0 or 1
    
    Returns: 
        DOFs of each node in array nodes.
    '''
    dofs = 2
    all_dofs = []
    mask = abs(nodes_dir[:, 1]) == 1
    all_dofs.extend(dofs * (nodes_dir[mask, 0] - 1))

    mask = abs(nodes_dir[:, 2]) == 1
    all_dofs.extend((dofs * nodes_dir[mask, 0]) - 1)

    all_dofs.sort()

    return np.array(all_dofs, dtype='int')

def get_load_vector(nelx, nely, force_matrix):
    """ Creates the force vector.

    Args:
        nelx (:obj:`int`): Number of elements on the x-axis.
        nely (:obj:`int`): Number of elements on the y-axis.
        force_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value.

    Returns:
        Loading vector.
    """
    ngl = 2 * ((nelx + 1) * (nely + 1))
    load_vector = np.zeros(ngl, dtype=complex)
    if len(force_matrix) > 1:
        force_matrix = force_matrix[np.argsort(force_matrix[:, 0])]       
    force_ind = get_dofs(force_matrix)
    load_vector[force_ind] = duplicate_force(force_matrix)

    return load_vector

def duplicate_force(force_matrix):
    """ Doubled force value for x direction = y direction = 1.

    Args:
        force_matrix (:obj:`numpy.array`): The columns are respectively node, x direction, y direction, force value.
    
    Returns:
        Force values.
    """
    mask = ((abs(force_matrix[:, 1]) == 1.) & (abs(force_matrix[:, 2]) == 1.)).nonzero()[0]
    force_values = force_matrix[:, 3]
    if mask.size != 0:
        force_values = np.insert(force_values, mask, force_matrix[mask, 3])
    # Change force sign
    aux = force_matrix[:, [1,2]].ravel()
    aux = aux[aux!=0]
    if (aux<0).any():
        force_values[aux<0] *= -1

    return force_values

def remove_dofs(nelx, nely, del_dofs):
    """ Delete specific DOFs from all DOFs.
    
    Args: 
        nelx (:obj:`int`): Number of elements on the x-axis.
        nely (:obj:`int`): Number of elements on the y-axis.
        del_dofs (:obj:`numpy.array`): Array with DOFs to be removed.
    
    Returns:
        Array without the DOFs removed.
    """
    dofs = np.arange((nelx+1) * (nely + 1) * 2)
    return np.delete(dofs, del_dofs)

def change_U_shape(disp_vector):
    ''' Transform displacement vector in matrix.
    
    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
    
    Returns: 
        Displacement of each axis.
    '''
    new_U = np.empty((int(len(disp_vector)/2), 2))
    new_U[:,0] = disp_vector[::2]
    new_U[:,1] = disp_vector[1::2] 
    return new_U

def apply_U(disp_vector, coord, factor):
    ''' Apply displacement to coordinates. 

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        coord (:obj:`numpy.array`): mesh coordinates.
        factor (:obj:`float`): Factor to deform the mesh.
    
    Returns: 
        Displaced mesh coordinates.
    '''
    new_coord = coord.copy()
    new_coord[:, 1:] += disp_vector * factor
    return new_coord

def freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, freq_rsp, delta, node_plot, load_vector, **kwargs):
    """ Get the displacement values for a specific node.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        ind_rows (:obj:`numpy.array`): Node indexes to make the assembly.
        ind_cols (:obj:`numpy.array`): Node indexes to make the assembly.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.  
        rho (:obj:`float`): Density.  
        alpha (:obj:`float`): Damping coefficient proportional to mass. 
        beta (:obj:`float`): Damping coefficient proportional to stiffness.  
        eta (:obj:`float`): Damping coefficient. 
        freq_rsp (:obj:`list`): Frequency range.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (:obj:`int`): Step between each calculation of the objective function. 
        node_plot (:obj:`int`): Node to salve the displacement.
        load_vector (:obj:`numpy.array`): Force.

    Returns:
        Displacement array.        
    """
    free_ind = None
    if kwargs.get('unrestricted_ind') is not None:
        free_ind = kwargs.get('unrestricted_ind')
    interval = np.arange(freq_rsp[0], freq_rsp[1] + 1, delta)
    vector_U = np.empty((len(interval)), dtype=complex)
    force_ind = get_dofs(node_plot)
    for n in range(len(interval)):
        disp_vector = solution2D(coord, connect, ind_rows, ind_cols, nelx, nely, E, v, rho, alpha, beta, eta, freq=interval[n], load_vector = load_vector, unrestricted_ind = free_ind)
        vector_U[n] = disp_vector[force_ind]
        print('It.', n)
        
    return vector_U