import gmsh 
import numpy as np    

def import_mesh(path):
    coord, connect = generate(path)
    return coord, connect

def generate(path):
    _initialize_gmsh(path)
    _set_gmsh_options()
    nodal_coordinates_matrix, connectivity_matrix = _get_matrices()
    _finalize_gmsh()
    return nodal_coordinates_matrix, connectivity_matrix

def _initialize_gmsh(path):
    gmsh.initialize('', False)
    gmsh.open(path)

def _set_gmsh_options():
    gmsh.option.setNumber('Mesh.Optimize', 1)
    gmsh.option.setNumber('Mesh.OptimizeNetgen', 0)
    gmsh.option.setNumber('Mesh.HighOrderOptimize', 0)
    gmsh.option.setNumber('Mesh.ElementOrder', 1)
    gmsh.option.setNumber('Mesh.Algorithm', 6)
    gmsh.option.setNumber('Mesh.RandomFactor', 1e-6)
    gmsh.option.setNumber('Geometry.Tolerance', 1e-4)

def _finalize_gmsh():
    gmsh.finalize()

def _get_matrices():
    # Algorithm 2
    for s in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(s[1])
        gmsh.model.mesh.setRecombine(s[0], s[1])
    gmsh.model.mesh.generate(2)        
    gmsh.model.mesh.removeDuplicateNodes()

    node_indexes, coords, _ = gmsh.model.mesh.getNodes(-1, -1, False)
    (unique, counts) = np.unique(node_indexes, return_counts=True)
    frequencies = np.asarray((unique, counts)).T
    _, element_indexes, connectivity = gmsh.model.mesh.getElements()
    
    map_nodes = dict(zip(node_indexes, np.arange(1, len(node_indexes)+1, 1)))
    nodal_coordinates_matrix = _get_nodal_coordinates_matrix(node_indexes, coords, map_nodes)
    
    nelx = len(nodal_coordinates_matrix[nodal_coordinates_matrix[:, 2] == nodal_coordinates_matrix[0, 2]]) - 1
    nely = len(nodal_coordinates_matrix[nodal_coordinates_matrix[:, 1] == nodal_coordinates_matrix[0, 1]]) - 1
    connectivity_matrix = _get_connectivity_matrix(element_indexes[1], connectivity[1], nelx, nely) #- (element_indexes[1][0] - 1) 
    return nodal_coordinates_matrix, connectivity_matrix
   
def _get_nodal_coordinates_matrix(indexes, coords, map_nodes):
    nodal_coordinates_matrix = np.zeros((len(indexes), 3), dtype=float)
    for i, (index, coord) in enumerate(zip(indexes, split_sequence(coords, 3))):
        x = mm_to_m(coord[0])
        y = mm_to_m(coord[1])
        nodal_coordinates_matrix[i,:] = [map_nodes[index], x, y]
    return nodal_coordinates_matrix

def _get_connectivity_matrix(total_el, connectivities, nelx, nely):
    indexes = np.arange(0, len(total_el))
    
    connectivity_matrix = np.zeros((len(total_el), 5), dtype=int)
    connectivity_matrix[:, 0] = indexes + 1
    connectivity_gmsh = split_sequence(connectivities, 4)
    
    indexes_gmsh = indexes.reshape(nely, nelx, order='F')
    indexes_fem = indexes.reshape(nely, nelx)
    for row_gmsh, row_fem in zip(indexes_gmsh, indexes_fem):
        for col_gmsh, col_fem in zip(row_gmsh, row_fem):
            connectivity_matrix[col_fem, 1:] = connectivity_gmsh[col_gmsh]
    return connectivity_matrix

def split_sequence(sequence, size):
    subsequences = []
    for start in range(0, len(sequence), size):
        end = start + size
        subsequence = sequence[start:end]
        subsequences.append(subsequence)
    return subsequences

def mm_to_m(m):
    return float(m) / 1000

def m_to_mm(m):
    return float(m) * 1000