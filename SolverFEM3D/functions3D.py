from time import time
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy import spatial

def shapeH8(rrx, ssx, ttx):
    """ Linear Shape Functions and Derivatives.

    Args:
        rrx (float): Local coordinate of the element on the Z-axis.
        ssx (float): Local coordinate of the element on the X-axis.
        ttx (float): Local coordinate of the element on the Y-axis.
    
    Returns:
        Tuple with the shape function and its derivative.
    """
    denominator = 8
    #shape functions 
    phi = np.zeros(8) 
    phi[0]=(1.-ssx)*(1.-ttx)*(1.-rrx)
    phi[1]=(1.+ssx)*(1.-ttx)*(1.-rrx)
    phi[2]=(1.+ssx)*(1.+ttx)*(1.-rrx)
    phi[3]=(1.-ssx)*(1.+ttx)*(1.-rrx)
    phi[4]=(1.-ssx)*(1.-ttx)*(1.+rrx)
    phi[5]=(1.+ssx)*(1.-ttx)*(1.+rrx)
    phi[6]=(1.+ssx)*(1.+ttx)*(1.+rrx)
    phi[7]=(1.-ssx)*(1.+ttx)*(1.+rrx)
    phi = phi/denominator
    #derivatives
    dphi = np.zeros((3,8))
    #
    dphi[0,0]=(1.-ssx)*(1.-ttx)*(-1.)
    dphi[0,1]=(1.+ssx)*(1.-ttx)*(-1.)
    dphi[0,2]=(1.+ssx)*(1.+ttx)*(-1.)
    dphi[0,3]=(1.-ssx)*(1.+ttx)*(-1.)
    dphi[0,4]=(1.-ssx)*(1.-ttx)*(1.)
    dphi[0,5]=(1.+ssx)*(1.-ttx)*(1.)
    dphi[0,6]=(1.+ssx)*(1.+ttx)*(1.)
    dphi[0,7]=(1.-ssx)*(1.+ttx)*(1.)
    #
    dphi[1,0]=(-1.)*(1.-ttx)*(1.-rrx)
    dphi[1,1]=(1.)*(1.-ttx)*(1.-rrx)
    dphi[1,2]=(1.)*(1.+ttx)*(1.-rrx)
    dphi[1,3]=(-1.)*(1.+ttx)*(1.-rrx)
    dphi[1,4]=(-1.)*(1.-ttx)*(1.+rrx)
    dphi[1,5]=(1.)*(1.-ttx)*(1.+rrx)
    dphi[1,6]=(1.)*(1.+ttx)*(1.+rrx)
    dphi[1,7]=(-1.)*(1.+ttx)*(1.+rrx)
    #
    dphi[2,0]=(1.-ssx)*(-1.)*(1.-rrx)
    dphi[2,1]=(1.+ssx)*(-1.)*(1.-rrx)
    dphi[2,2]=(1.+ssx)*(1.)*(1.-rrx)
    dphi[2,3]=(1.-ssx)*(1.)*(1.-rrx)
    dphi[2,4]=(1.-ssx)*(-1.)*(1.+rrx)
    dphi[2,5]=(1.+ssx)*(-1.)*(1.+rrx)
    dphi[2,6]=(1.+ssx)*(1.)*(1.+rrx)
    dphi[2,7]=(1.-ssx)*(1.)*(1.+rrx)
    #
    dphi = dphi/denominator
    #
    return phi, dphi
      
def matricesH8(ee, coord, connect, E, v, rho):
    """ H8 stiffness and mass matrices.

    Args:
        ee (int): Element.
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.
        rho (float): Density.

    Returns:
        Tuple with stiffness and mass matrices.
    """
    # constitutive info
    CTTV = np.zeros((6,6))
    tempc=E/((1.+v)*(1.-2.*v))
    tempn=(1.-2.*v)/2.
    tempt=1.-v
    CTTV[0,0]=tempt
    CTTV[0,1]=v
    CTTV[0,2]=v
    CTTV[1,0]=v
    CTTV[1,1]=tempt
    CTTV[1,2]=v
    CTTV[2,0]=v
    CTTV[2,1]=v
    CTTV[2,2]=tempt
    CTTV[3,3]=tempn
    CTTV[4,4]=tempn
    CTTV[5,5]=tempn
    CTTV=tempc*CTTV
    # integration points
    nint, con, wps = 8, 1/np.sqrt(3), 1
    pint = np.zeros((8,3)) + con
    pint[0,0]=-con
    pint[0,1]=-con
    pint[0,2]=-con
    pint[1,1]=-con
    pint[1,2]=-con
    pint[2,2]=-con
    pint[3,0]=-con
    pint[3,2]=-con
    pint[4,0]=-con
    pint[4,1]=-con
    pint[5,1]=-con
    pint[7,0]=-con
    # preallocating elementary matrices
    Ke, Me = 0, 0
    AUJJ = np.zeros((3,3))
    B = np.zeros((6,24))
    N = np.zeros((3,24))
    # integration
    for i in range(nint):
        rrx, ssx, ttx = pint[i, 0], pint[i, 1], pint[i, 2]
        phi, dphi = shapeH8(ssx,ttx,rrx)
        ie = connect[ee,1:]-1
        dxdydz = dphi@coord[ie, 1:4] 
        # note: dxds, dyds, dzds, dxdt, dydt, dzdt, dxdr, dydr, dzdr 
        # = dxdydz[0,0], dxdydz[0,1], dxdydz[0,2], dxdydz[1,0],dxdydz[1,1],dxdydz[1,2], dxdydz[2,0],dxdydz[2,1],dxdydz[2,2]  
        JAC = np.array([[dxdydz[0,0], dxdydz[0,1], dxdydz[0,2]],
                        [dxdydz[1,0], dxdydz[1,1], dxdydz[1,2]],
                        [dxdydz[2,0], dxdydz[2,1], dxdydz[2,2]]], dtype=float) # JAC = np.array([[dxds, dyds],[dxdt, dydt]], dtype=float)
        detJAC = (JAC[0,0] * JAC[1,1] * JAC[2,2] + 
                JAC[0,1] * JAC[1,2] * JAC[2,0] + 
                JAC[0,2] * JAC[1,0] * JAC[2,1]) - \
                (JAC[2,0] * JAC[1,1] * JAC[0,2] + 
                JAC[2,1] * JAC[1,2] * JAC[0,0] + 
                JAC[2,2] * JAC[1,0] * JAC[0,1])
        # adj(JAC)
        AUJJ[0,0]= 1 * ((JAC[1,1] * JAC[2,2]) - (JAC[2,1] * JAC[1,2]))
        AUJJ[1,0]= -1 * ((JAC[1,0] * JAC[2,2]) - (JAC[1,2] * JAC[2,0]))
        AUJJ[2,0]= 1 * ((JAC[1,0] * JAC[2,1]) - (JAC[1,1] * JAC[2,0]))
        AUJJ[0,1]= -1 * ((JAC[0,1] * JAC[2,2]) - (JAC[0,2] * JAC[2,1]))
        AUJJ[1,1]= 1 * ((JAC[0,0] * JAC[2,2]) - (JAC[0,2] * JAC[2,0]))
        AUJJ[2,1]= -1 * ((JAC[0,0] * JAC[2,1]) - (JAC[0,1] * JAC[2,0]))
        AUJJ[0,2]= 1 * ((JAC[0,1] * JAC[1,2]) - (JAC[0,2] * JAC[1,1]))
        AUJJ[1,2]= -1 * ((JAC[0,0] * JAC[1,2]) - (JAC[0,2] * JAC[1,0]))
        AUJJ[2,2]= 1 * ((JAC[0,0] * JAC[1,1]) - (JAC[0,1] * JAC[1,0]))
        #Inverse Jacobian
        iJAC = (1/detJAC) * AUJJ # np.linalg.inv(JAC) 
        #
        dphi_t = iJAC @ dphi
        #
        for iii in range(8):
            B[0,3*(iii)+0]=dphi_t[0,iii]
            B[0,3*(iii)+1]=0.
            B[0,3*(iii)+2]=0.
            B[1,3*(iii)+0]=0.
            B[1,3*(iii)+1]=dphi_t[1,iii]
            B[1,3*(iii)+2]=0.
            B[2,3*(iii)+0]=0.
            B[2,3*(iii)+1]=0.
            B[2,3*(iii)+2]=dphi_t[2,iii]
            B[3,3*(iii)+0]=dphi_t[1,iii]
            B[3,3*(iii)+1]=dphi_t[0,iii]
            B[3,3*(iii)+2]=0.              
            B[4,3*(iii)+0]=0.
            B[4,3*(iii)+1]=dphi_t[2,iii]
            B[4,3*(iii)+2]=dphi_t[1,iii]
            B[5,3*(iii)+0]=dphi_t[2,iii]
            B[5,3*(iii)+1]=0.
            B[5,3*(iii)+2]=dphi_t[0,iii]
        #
        for iii in range(8):
            N[0,3*(iii)+0]=phi[iii]
            N[0,3*(iii)+1]=0.
            N[0,3*(iii)+2]=0.
            N[1,3*(iii)+0]=0.
            N[1,3*(iii)+1]=phi[iii]
            N[1,3*(iii)+2]=0.
            N[2,3*(iii)+0]=0.
            N[2,3*(iii)+1]=0.
            N[2,3*(iii)+2]=phi[iii]
      
        #
        Ke += B.T@(CTTV@B)*(detJAC*wps)
        Me += rho*N.T@N*(detJAC*wps)
        #
    return Ke, Me

def regularmeshH8(nelx, nely, nelz, lx, ly, lz):
    """ Create a regular H8 mesh.

    Args:
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        nelz (int): Number of elements on the Z-axis.
        lx (int): X-axis length.
        ly (int): Y-axis length.
        lz (int): Z-axis length.
    
    Returns:
        Tuple with the coordinate matrix, connectivity, and the indexes of each node.    
    """
    x, y, z = np.linspace(0, lx, num=nelx + 1), np.linspace(0, ly, num=nely + 1), np.linspace(0, lz, num=nelz + 1)
    nx, ny, nz = len(x), len(y), len(z)
    mat_x = (x.reshape(nx, 1)@np.ones((1, ny*nz))).T
    mat_y = y.reshape(ny, 1)@np.ones((1, nx))
    mat_z = z.reshape(nz, 1)@np.ones((1, nx*ny))
    x_t, y_t, z_t = mat_x.flatten(), np.tile(mat_y.flatten(), nz), mat_z.flatten()
    ind_coord = np.arange(1, (nz)* nx * ny + 1, 1, dtype=int)
    coord = (np.array([ind_coord, x_t, y_t, z_t])).T
    # processing of connectivity matrix
    ind_connect = np.arange(1, nelz * nelx * nely + 1, dtype=int)
    mat_aux = ind_connect.reshape(nely, nelx, nelz)
    a = np.arange(0, nely * nelz, 1)
    for ind in range(nely, len(a), nely):
        a[ind:] += nelx + 1
    c = (a.reshape(len(a),1)@np.ones((1, nelx))).reshape(nely, nelx, nelz)
    b = (mat_aux + c).flatten()

    connect = np.array([ind_connect, b+(nelx+1), b, b+1, b+(nelx+2), \
                        b+(nelx+1)*(nely+1)+(nelx+1), b+(nelx+1)*(nely+1), \
                        b+1+(nelx+1)*(nely+1), b+(nelx+1)*(nely+1)+(nelx+2)], dtype=int).T
    # processing the dofs indices (rows and columns) for assembly
    dofs, edofs = 3, 24
    ind_dofs = (np.array([dofs*connect[:,1]-1, dofs*connect[:,1], dofs*connect[:,1]+1,
                        dofs*connect[:,2]-1, dofs*connect[:,2], dofs*connect[:,2]+1,
                        dofs*connect[:,3]-1, dofs*connect[:,3], dofs*connect[:,3]+1,
                        dofs*connect[:,4]-1, dofs*connect[:,4], dofs*connect[:,4]+1,
                        dofs*connect[:,5]-1, dofs*connect[:,5], dofs*connect[:,5]+1,
                        dofs*connect[:,6]-1, dofs*connect[:,6], dofs*connect[:,6]+1,
                        dofs*connect[:,7]-1, dofs*connect[:,7], dofs*connect[:,7]+1,
                        dofs*connect[:,8]-1, dofs*connect[:,8], dofs*connect[:,8]+1], dtype=int)-2).T
    
    vect_indices = ind_dofs.flatten()
    ind_rows = ((np.tile(vect_indices, (edofs,1))).T).flatten()
    ind_cols = (np.tile(ind_dofs, edofs)).flatten()
    
    return coord, connect, ind_rows, ind_cols
   
def solution3D(coord, connect, ind_rows, ind_cols, nelx, nely, nelz, E, v, rho, alpha, beta, eta, freq, timing=False, **kwargs):
    """ Assembly and solution.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        ind_rows (numpy.array): Node indexes to make the assembly.
        ind_cols (numpy.array): Node indexes to make the assembly.
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        nelz (int): Number of elements on the Z-axis.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.  
        rho (float): Density.  
        alpha (float): Damping coefficient proportional to mass. 
        beta (float): Damping coefficient proportional to stiffness.  
        eta (float): Damping coefficient. 
        freq (int): Analyzed frequency.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.

    Returns:
        Displacement array.
    """
    #
    t01 = time()
    ngl = 3 * ((nelx + 1) * (nely + 1) * (nelz + 1))
    data_k = np.zeros((nelx * nely * nelz, 576), dtype=float)
    data_m = np.zeros((nelx * nely * nelz, 576), dtype=float)
    #
    for el in range(nelx * nely * nelz):
        Ke, Me = matricesH8(el, coord, connect, E, v, rho)
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
    w = 2 * np.pi * freq
    if w == 0:
        w = 1e-12
    damp_matrix = alpha * mass_matrix + (beta + eta/w)*stif_matrix
    if kwargs.get('load_vector') is not None:
        load_vector = kwargs.get('load_vector')  
    #
    if kwargs.get('unrestricted_ind') is not None:
        free_ind = kwargs.get('unrestricted_ind')
        F = load_vector[free_ind]
        K = stif_matrix[free_ind, :][:, free_ind]
        M = mass_matrix[free_ind, :][:, free_ind]
        C = damp_matrix[free_ind, :][:, free_ind]
        #Kd = -(w**2)*M + K  
        Kd = -(w**2)*M + 1j*w*C + K
        U = spsolve(Kd,F)
        disp_vector = np.zeros((ngl), dtype=complex)
        disp_vector[free_ind] = U
    else:
        #Kd = -(w**2)*mass_matrix + stif_matrix
        Kd = -(w**2)*mass_matrix + 1j*w*damp_matrix + stif_matrix
        disp_vector = spsolve(Kd, load_vector)
    tf2 = time()
    if timing:
        print("Time to solve the harmonic analysis problem: " + str(round((tf2 - t02),6)) + '[s]')
        print("Total time elapsed in solution: " + str(round((tf2 - t01),6)) + '[s]')
    #
    return disp_vector

def freqresponse(coord, connect, ind_rows, ind_cols, nelx, nely, nelz, E, v, rho, alpha, beta, eta, freq_range, delta, node_plot, load_vector, **kwargs):
    """ Get the displacement values for a specific node.

    Args:
        coord (numpy.array): Coordinates of the element.
        connect (numpy.array): Element connectivity.
        ind_rows (numpy.array): Node indexes to make the assembly.
        ind_cols (numpy.array): Node indexes to make the assembly.
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        nelz (int): Number of elements on the Z-axis.
        E (float): Elastic modulus.
        v (float): Poisson's ratio.  
        rho (float): Density.  
        alpha (float): Damping coefficient proportional to mass. 
        beta (float): Damping coefficient proportional to stiffness.  
        eta (float): Damping coefficient. 
        freq_rsp (list): Frequency range.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (int): Step between each calculation of the objective function. 
        node_plot (int): Node to salve the displacement.
        load_vector (numpy.array): Force.

    Returns:
        Displacement array.        
    """
    free_ind = None
    if kwargs.get('unrestricted_ind') is not None:
        free_ind = kwargs.get('unrestricted_ind')
    interval = np.arange(freq_range[0], freq_range[1] + 1, delta)
    vector_U = np.empty((len(interval)), dtype=complex)
    force_ind = get_dofs(node_plot)
    for n in range(len(interval)):
        disp_vector = solution3D(coord, connect, ind_rows, ind_cols, nelx, nely, nelz, E, v, rho, alpha, beta, eta, freq=interval[n], load_vector = load_vector, unrestricted_ind = free_ind)
        vector_U[n] = disp_vector[force_ind]
        print('It.', n)
    return vector_U

def change_U_shape(disp_vector):
    ''' Transform displacement vector in matrix.
    
    Args:
        disp_vector (numpy.array): Displacement.
    
    Returns: 
        Displacement of each axis.
    '''
    new_U = np.empty((int(len(disp_vector)/3), 3))
    new_U[:,0] = disp_vector[0:len(disp_vector):3]
    new_U[:,1] = disp_vector[1:len(disp_vector):3]
    new_U[:,2] = disp_vector[2:len(disp_vector):3]

    return new_U

def apply_U(disp_vector, coord, factor):
    ''' Apply displacement to coordinates. 

    Args:
        disp_vector (numpy.array): Displacement.
        coord (numpy.array): mesh coordinates.
        factor (float): Factor to deform the mesh.
    
    Returns: 
        Displaced mesh coordinates.
    '''
    new_coord = coord.copy()
    new_coord[:, 1:] += disp_vector * factor

    return new_coord

def get_nodes_by_coord(coord, coord_user):
    ''' Get node numbers by coordinate.

    Args:
        coord (numpy.array): mesh coordinates.
        coord_user (numpy.array): user coordinates.
    
    Returns: 
        Nodes of the coordinates provided.
    '''
    mytree = spatial.cKDTree(coord[:, [1,2,3]])
    _, ind_nodes = mytree.query(coord_user)
    nodes = coord[ind_nodes, 0]
    return nodes

def get_nodes1d(coord, coord_user, eps, column):
    ''' Get node numbers that are equal to coord.

    Args:
        coord (numpy.array): mesh coordinates.
        coord_user (numpy.array): coordinates in one direction (x, y or z).
        eps (float): Acceptable margin of difference.
        column (int): Direction to compare (x, y or z).

    Returns:
        Nodes.
    '''
    dif = np.abs(coord[:, column] - coord_user)
    mask = dif < eps
    return (coord[mask, 0]).astype('int')

def get_dofs(nodes_direct):
    ''' Get DOFs that meet the specified direction.

    Args:
        nodes_dir (numpy.array): [nodes numbers, x_direction, y_direction, z_direction].
            x_direction, y_direction and z_direction can be -1, 0 or 1.
    
    Returns: 
        DOFs of each node in array nodes.
    '''
    dofs = 3
    all_dofs = []
    mask = abs(nodes_direct[:, 1]) == 1
    all_dofs.extend((dofs * nodes_direct[mask, 0]) - 3)

    mask = abs(nodes_direct[:, 2]) == 1
    all_dofs.extend((dofs * nodes_direct[mask, 0]) - 2)

    mask = abs(nodes_direct[:, 3]) == 1
    all_dofs.extend((dofs * nodes_direct[mask, 0]) - 1)

    all_dofs.sort()

    return np.array(all_dofs, dtype='int')

def remove_dofs(nelx, nely, nelz, del_dofs):
    """ Delete specific DOFs from all DOFs.
    
    Args: 
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        nelz (int): Number of elements on the Z-axis.
        del_dofs (numpy.array): Array with DOFs to be removed.
    
    Returns:
        Array without the DOFs removed.
    """
    dofs = np.arange((nelx+1) * (nely + 1) * nelz * 3)
    return np.delete(dofs, del_dofs)

def duplicate_force(load_matrix):
    """ Doubled force value.

    Args:
        load_matrix (numpy.array): Force.
    
    Returns:
        Load values.
    """
    load_vector = load_matrix[:, -1]
    aux = np.sum(abs(load_matrix[:, 1:-1]), axis=1)
    idx_dup = np.where((aux>1) & (aux<3))[0]
    idx_trip = np.repeat(np.where(aux>2), 2)
    if idx_dup.size != 0 or idx_trip.size != 0:
        load_vector = np.insert(load_vector, np.hstack([idx_dup, idx_trip]), np.hstack([load_matrix[idx_dup, -1], load_matrix[idx_trip, -1]]))
    # Change force sign
    aux = load_matrix[:, [1, 2, 3]].ravel()
    aux = aux[aux!=0]
    if (aux<0).any():
        load_vector[aux < 0] *= -1

    return load_vector

def get_load_vector(nelx, nely, nelz, load_matrix):
    """ Creates the load vector.

    Args:
        nelx (int): Number of elements on the X-axis.
        nely (int): Number of elements on the Y-axis.
        nelz (int): Number of elements on the Z-axis.
        load_matrix (numpy.array): Force.

    Returns:
        Loading vector.
    """
    ngl = 3 * ((nelx + 1) * (nely + 1) * (nelz + 1))
    load_vector = np.zeros(ngl, dtype=complex)
    if len(load_matrix) > 1:
        load_matrix = load_matrix[np.argsort(load_matrix[:, 0])]  
    force_ind = get_dofs(load_matrix)
    load_vector[force_ind] = duplicate_force(load_matrix)

    return load_vector

def get_max_min_freq(freq_range, delta, disp_vector):
    """ Get the frequency with the minimum and maximum displacement, respectively.

    Args:
        freq_range (list): The initial and final frequency.
        delta (int): step between each calculation of the displacement.
        disp_vector (numpy.array): Displacement.
    
    Returns:
        A tuple with the minimum and maximum frequency.
    """
    x = np.arange(freq_range[0], freq_range[1] + 1, delta)
    y = 10 * np.log10(abs(disp_vector))
    return (x[np.where(y == y.min())], x[np.where(y == y.max())])