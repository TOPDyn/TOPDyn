from time import time
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy import spatial

def shapeH8(rrx, ssx, ttx):
    """ Linear Shape Functions and Derivatives.

    Args:
        rrx (:obj:`float`): Local coordinate of the element on the Z-axis.
        ssx (:obj:`float`): Local coordinate of the element on the X-axis.
        ttx (:obj:`float`): Local coordinate of the element on the Y-axis.
    
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
    
    dphi[0,0]=(1.-ssx)*(1.-ttx)*(-1.)
    dphi[0,1]=(1.+ssx)*(1.-ttx)*(-1.)
    dphi[0,2]=(1.+ssx)*(1.+ttx)*(-1.)
    dphi[0,3]=(1.-ssx)*(1.+ttx)*(-1.)
    dphi[0,4]=(1.-ssx)*(1.-ttx)*(1.)
    dphi[0,5]=(1.+ssx)*(1.-ttx)*(1.)
    dphi[0,6]=(1.+ssx)*(1.+ttx)*(1.)
    dphi[0,7]=(1.-ssx)*(1.+ttx)*(1.)
    
    dphi[1,0]=(-1.)*(1.-ttx)*(1.-rrx)
    dphi[1,1]=(1.)*(1.-ttx)*(1.-rrx)
    dphi[1,2]=(1.)*(1.+ttx)*(1.-rrx)
    dphi[1,3]=(-1.)*(1.+ttx)*(1.-rrx)
    dphi[1,4]=(-1.)*(1.-ttx)*(1.+rrx)
    dphi[1,5]=(1.)*(1.-ttx)*(1.+rrx)
    dphi[1,6]=(1.)*(1.+ttx)*(1.+rrx)
    dphi[1,7]=(-1.)*(1.+ttx)*(1.+rrx)
    
    dphi[2,0]=(1.-ssx)*(-1.)*(1.-rrx)
    dphi[2,1]=(1.+ssx)*(-1.)*(1.-rrx)
    dphi[2,2]=(1.+ssx)*(1.)*(1.-rrx)
    dphi[2,3]=(1.-ssx)*(1.)*(1.-rrx)
    dphi[2,4]=(1.-ssx)*(-1.)*(1.+rrx)
    dphi[2,5]=(1.+ssx)*(-1.)*(1.+rrx)
    dphi[2,6]=(1.+ssx)*(1.)*(1.+rrx)
    dphi[2,7]=(1.-ssx)*(1.)*(1.+rrx)
    
    dphi = dphi/denominator
    return phi, dphi
      
def matricesH8(ee, coord, connect, E, v, rho):
    """ H8 stiffness and mass matrices.

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
        
        dphi_t = iJAC @ dphi
        
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
      
        Ke += B.T@(CTTV@B)*(detJAC*wps)
        Me += rho*N.T@N*(detJAC*wps)
    return Ke, Me

def regularmeshH8(nelx, nely, nelz, lx, ly, lz):
    """ Creates a regular H8 mesh.

    Args:
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        nelz (:obj:`int`): Number of elements on the Z-axis.
        lx (:obj:`float`): X-axis length.
        ly (:obj:`float`): Y-axis length.
        lz (:obj:`float`): Z-axis length.
    
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
    return coord, connect

def generate_ind_rows_cols(connect):
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
    return ind_rows, ind_cols

def turn_into_np(dicti_matrix, force):
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
                aux = [dict_row["x_coord"], dict_row["y_coord"], dict_row["z_coord"], dict_row["x_direc"], dict_row["y_direc"], dict_row["z_direc"], dict_row["force"]]  
            else:
                aux = [dict_row["x_coord"], dict_row["y_coord"], dict_row["z_coord"], dict_row["constrain_disp_x"], dict_row["constrain_disp_y"], dict_row["constrain_disp_z"]] 
            index_by_coord.append(i)
        
        elif "eps" in dict_row:
            if force:
                aux = [dict_row["coord"], dict_row["axis"], dict_row["eps"], dict_row["x_direc"], dict_row["y_direc"], dict_row["z_direc"], dict_row["force"]] 
            else:
                aux = [dict_row["coord"], dict_row["axis"], dict_row["eps"], dict_row["constrain_disp_x"], dict_row["constrain_disp_y"], dict_row["constrain_disp_z"]]
            index_by_col.append(i)
        
        matrix.append(aux)

    matrix_F = np.array(matrix)
    return index_by_coord, index_by_col, matrix_F

def get_matrices(matrix, coord, force):
    """ Gets the force matrix or the contraint matrix.

    Args:
        matrix (:obj:`list`): List passed by the user. 
        coord (:obj:`numpy.array`): Coordinates of the element.
        force (:obj:`bool`): True if encountering force matrix.
    Returns:
        force_matrix or restri_matrix.
    """
    if force:
        index_by_coord, index_by_col, np_matrix = turn_into_np(matrix, force=True)
        nodes_coord, nodes_col = get_nodes(coord, np_matrix, index_by_coord, index_by_col)
        nodes_matrix = get_matrix(nodes_coord, nodes_col, np_matrix, index_by_coord, index_by_col, [1,2,3,4], [3,4,5,6], 5)
    else:
        index_by_coord, index_by_col, np_matrix = turn_into_np(matrix, force=False)
        nodes_coord, nodes_col = get_nodes(coord, np_matrix, index_by_coord, index_by_col)
        nodes_matrix = get_matrix(nodes_coord, nodes_col, np_matrix, index_by_coord, index_by_col, [1,2,3], [3,4,5], 4)
    return nodes_matrix.astype(int)

def get_nodes(coord, np_matrix, index_by_coord, index_by_col):
    """ Gets nodes by a coordinate or a column.

    Args:
        coord (:obj:`numpy.array`): Coordinates of the element.
        np_matrix (:obj:`numpy.array`): List passed to an array.
        index_by_coord (:obj:`list`): indices of elements passed by coordinates.
        index_by_col (:obj:`list`): indices of elements passed by columns.
        ind (:obj:`int`): The column that has the margin of error.
    Returns:
        Nodes.
    """
    nodes_coord = []
    nodes_col = []

    if len(index_by_coord) > 0:   
        nodes_coord = get_nodes_by_coord(coord, np_matrix[np.ix_(index_by_coord, [0,1,2])])

    if len(index_by_col) > 0:
        for index in index_by_col:
            aux = get_nodes1d(coord, np_matrix[index, 0], np_matrix[index, 2], int(np_matrix[index, 1]))
            nodes_col.append(aux)
    return nodes_coord, nodes_col 

def get_matrix(nodes_coord, nodes_col, np_matrix, index_by_coord, index_by_col, ind1, ind2, total_cols):
    """ Creates the node matrix.

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

def stif_mass_matrices(nelx, nely, nelz, coord, connect, ind_rows, ind_cols, E, v, rho, timing):
    """ Calculates global matrices.

    Args:
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        nelz (:obj:`int`): Number of elements on the Z-axis.
        coord (:obj:`numpy.array`): Coordinates of the element.
        connect (:obj:`numpy.array`): Element connectivity.
        ind_rows (:obj:`numpy.array`): Node indexes to make the assembly.
        ind_cols (:obj:`numpy.array`): Node indexes to make the assembly.
        E (:obj:`float`): Elastic modulus.
        v (:obj:`float`): Poisson's ratio.  
        rho (:obj:`float`): Density.  
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.

    Returns:
        Global matrices
    """
    tm = time()
    ngl = 3 * ((nelx + 1) * (nely + 1) * (nelz + 1))
    data_k = np.zeros((nelx * nely * nelz, 576), dtype=float)
    data_m = np.zeros((nelx * nely * nelz, 576), dtype=float)
    
    for el in range(nelx * nely * nelz):
        Ke, Me = matricesH8(el, coord, connect, E, v, rho)
        data_k[el,:] = Ke.flatten() 
        data_m[el,:] = Me.flatten() 
    
    data_k = data_k.flatten()
    data_m = data_m.flatten()
    stif_matrix = csc_matrix((data_k, (ind_rows, ind_cols)), shape=(ngl, ngl))
    mass_matrix = csc_matrix((data_m, (ind_rows, ind_cols)), shape=(ngl, ngl))
    tmf = time()
    if timing:
        print("Time to assembly global matrices: " + str(round((tmf - tm), 6)) + '[s]')
    return stif_matrix, mass_matrix

def get_damp_matrix(mass_matrix, stif_matrix, freq, alpha, beta, eta):
    """ Calculates damping matrix.

    Args:
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        freq (:obj:`int`): Analyzed frequency.
        alpha (:obj:`float`): Damping coefficient proportional to mass. 
        beta (:obj:`float`): Damping coefficient proportional to stiffness.  
        eta (:obj:`float`): Damping coefficient. 

    Returns:
        Damping matrix. 
    """
    w = 2 * np.pi * freq
    if w == 0:
        w = 1e-12
    return alpha * mass_matrix + (beta + eta/w)*stif_matrix

def get_displacement(load_vector, free_ind, stif_matrix, mass_matrix, damp_matrix, freq, ngl, timing):
    """ Harmonic solution.

    Args:
        load_vector (:obj:`numpy.array`): Load.
        free_ind (:obj:`numpy.array`): Free dofs. 
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        damp_matrix (:obj:`numpy.array`): Damping matrix.
        freq (:obj:`int`): Analyzed frequency.
        ngl (:obj:`int`): Degrees of freedom.
        timing (:obj:`bool`, optional): If True shows the process optimization time. Defaults to False.

    Returns:
        Displacement vector.
    """  
    tf = time()
    w = 2 * np.pi * freq
    if free_ind is not None:
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
        print("Time to solve the harmonic analysis problem: " + str(round((tf2 - tf), 6)) + '[s]')
    return disp_vector

def freqresponse(stif_matrix, mass_matrix, load_vector, free_ind, ngl, alpha, beta, eta, freq_range, node_plot):
    """ Frequency response.

    Args:
        stif_matrix (:obj:`numpy.array`): Stiffness matrix.
        mass_matrix (:obj:`numpy.array`): Mass matrix.
        load_vector (:obj:`numpy.array`): Load.
        free_ind (:obj:`numpy.array`): Free dofs. 
        ngl (:obj:`int`): Degrees of freedom.             
        alpha (:obj:`float`): Damping coefficient proportional to mass. 
        beta (:obj:`float`): Damping coefficient proportional to stiffness.  
        eta (:obj:`float`): Damping coefficient. 
        freq_range (:obj:`list`): Frequency range.
            
            * First value is the minimum frequency.
            * Second value is the maximum frequency.
            * Third value is the step between each calculation of the objective function. 
        node_plot (:obj:`int`): Node to calculates the displacement.

    Returns:
        Displacement.
    """
    interval = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
    
    vector_U = np.empty((len(interval)), dtype=complex)
    force_ind = get_dofs(node_plot)
    
    for n, f in enumerate(interval):#range(len(interval)):
        damp_matrix = get_damp_matrix(mass_matrix, stif_matrix, f, alpha, beta, eta)
        disp_vector = get_displacement(load_vector, free_ind, stif_matrix, mass_matrix, damp_matrix, f, ngl, timing=False)      
        vector_U[n] = disp_vector[force_ind]
        print('It.', n)
    np.savetxt("U_errado", vector_U)
    return vector_U

def change_U_shape(disp_vector):
    """ Transforms displacement vector in matrix.
    
    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
    
    Returns: 
        Displacement of each axis.
    """
    new_U = np.empty((int(len(disp_vector)/3), 3))
    new_U[:,0] = disp_vector[0:len(disp_vector):3]
    new_U[:,1] = disp_vector[1:len(disp_vector):3]
    new_U[:,2] = disp_vector[2:len(disp_vector):3]

    return new_U

def apply_U(disp_vector, coord, factor):
    """ Applies displacement to coordinates. 

    Args:
        disp_vector (:obj:`numpy.array`): Displacement.
        coord (:obj:`numpy.array`): mesh coordinates.
        factor (:obj:`float`): Factor to deform the mesh.
    
    Returns: 
        Displaced mesh coordinates.
    """
    new_coord = coord.copy()
    new_coord[:, 1:] += disp_vector * factor

    return new_coord

def get_nodes_by_coord(coord, coord_user):
    """ Gets node numbers by coordinate.

    Args:
        coord (:obj:`numpy.array`): mesh coordinates.
        coord_user (:obj:`numpy.array`): user coordinates.
    
    Returns: 
        Nodes of the coordinates provided.
    """
    mytree = spatial.cKDTree(coord[:, [1,2,3]])
    _, ind_nodes = mytree.query(coord_user)
    nodes = coord[ind_nodes, 0]
    return nodes

def get_nodes1d(coord, coord_user, eps, column):
    """ Gets node numbers that are equal to coord.

    Args:
        coord (:obj:`numpy.array`): mesh coordinates.
        coord_user (:obj:`numpy.array`): coordinates in one direction (x, y or z).
        eps (:obj:`float`): Acceptable margin of difference.
        column (:obj:`int`): Direction to compare (x, y or z).

    Returns:
        Nodes.
    """
    dif = np.abs(coord[:, column] - coord_user)
    mask = dif < eps
    return (coord[mask, 0]).astype('int')

def get_dofs(nodes_direct):
    """ Gets DOFs that meet the specified direction.

    Args:
        nodes_dir (:obj:`numpy.array`): [nodes numbers, x_direction, y_direction, z_direction].
            x_direction, y_direction and z_direction can be -1, 0 or 1.
    
    Returns: 
        DOFs of each node in array nodes.
    """
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
    """ Deletes specific DOFs from all DOFs.
    
    Args: 
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        nelz (:obj:`int`): Number of elements on the Z-axis.
        del_dofs (:obj:`numpy.array`): Array with DOFs to be removed.
    
    Returns:
        Array without the DOFs removed.
    """
    dofs = np.arange((nelx+1) * (nely + 1) * nelz * 3)
    return np.delete(dofs, del_dofs)

def duplicate_force(load_matrix):
    """ Doubled force value.

    Args:
        load_matrix (:obj:`numpy.array`): Load.
    
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

def get_load_vector(ngl, load_matrix):
    """ Creates the load vector.

    Args:
        ngl (:obj:`int`): Degrees of freedom.
        load_matrix (:obj:`numpy.array`): Load.

    Returns:
        Loading vector.
    """
    load_vector = np.zeros(ngl, dtype=complex)
    if len(load_matrix) > 1:
        load_matrix = load_matrix[np.argsort(load_matrix[:, 0])]  
    force_ind = get_dofs(load_matrix)
    load_vector[force_ind] = duplicate_force(load_matrix)
    return load_vector

def get_max_min_freq(freq_range, delta, disp_vector):
    """ Gets the frequency with the minimum and maximum displacement, respectively.

    Args:
        freq_range (:obj:`list`): The initial and final frequency.
        delta (:obj:`int`): step between each calculation of the displacement.
        disp_vector (:obj:`numpy.array`): Displacement.
    
    Returns:
        A tuple with the minimum and maximum frequency.
    """
    x = np.arange(freq_range[0], freq_range[1] + 1, delta)
    y = 10 * np.log10(abs(disp_vector))
    return (x[np.where(y == y.min())], x[np.where(y == y.max())])