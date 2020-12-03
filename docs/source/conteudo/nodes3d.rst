Constrain nodes displacements in 3D
================


Figure 10 shows a mesh with *nelx = 2*, *nely = 2*, *nelz = 1*, and how the nodes are numbered. 


.. figure:: /conteudo/images/3d/mesh3d.JPG
   :scale: 50 %
   :align: center

   Figure 10: Mesh 3D with numbered nodes.

Constrain X, Y and Z displacements of the nodes in the plane X = 0
--------------------

A mask is created that returns an array with False and True values to select the nodes whose X coordinates in the coordinate matrix (column = 1) are equal to zero.

Then a matrix with four columns is created with the constrained nodes. The first column is the node number, the second column designates whether the node is constrained on X displacement (*1: True*, *0: False*), the third designates whether the node is constrained on Y displacement (*1: True*, *0: False*), and the fourth designates whether the node is constrained on Z displacement (*1: True*, *0: False*). 

In this example, the second, third, and fourth columns are equal to 1. The number of rows is defined seeing that the objective is to constrain all nodes of a rectangle face with a fixed X value. It will be :math:`(nely + 1) * (nelz + 1)`.

.. code-block:: python

    mask = coord[:, 1] == 0.
    restri_matrix = np.empty(((nely + 1) * (nelz + 1), 4))
    restri_matrix[:, 0] = (coord[mask, 0]).astype('int')
    restri_matrix[:, [1,2,3]] = np.ones(((nely + 1) * (nelz + 1), 3), dtype='int')

The resulting mesh is shown in figure 11.

.. figure:: /conteudo/images/3d/ex1.PNG
   :scale: 50 %
   :align: center

   Figure 11: Mesh with X, Y and Z axes constrained of the nodes in the plane X=0.

Constrain X displacement of the nodes in the plane Z = 0.3
--------------------

In this example, the second column is equal to 1, and the third and fourth columns are equal to 0. The number of rows is defined seeing that the objective is to constrain all nodes of a rectangle face with a fixed Z value. It will be :math:`(nely + 1) * (nelz + 1)`.

.. code-block:: python

    mask = coord[:, 3] == 0.3
    restri_matrix = np.empty(((nelx + 1) * (nely + 1), 4))
    restri_matrix[:, 0] = (coord[mask, 0]).astype('int')
    restri_matrix[:, 1] = np.ones(((nelx + 1) * (nely + 1)), dtype='int')
    restri_matrix[:, [2,3]] = np.zeros(((nelx + 1) * (nely + 1), 2), dtype='int')


The resulting mesh is shown in figure 12.

.. figure:: /conteudo/images/3d/ex2.PNG
   :scale: 50 %
   :align: center

   Figure 12: Mesh with X-axis constrained of the nodes in the plane Z=0.3.

Constrain X, Y and Z displacements of the nodes at coordinates: A = (1, 0.5, 0.3) and B = (1, 0.5, 0)
--------------------

Function :meth:`SolverFEM3D.functions3D.get_nodes_by_coord` works like 2D, but in 3D the Z-axis is added. The second, third, and fourth columns of the constrained nodes matrix are equal to 1.

.. code-block:: python

    restri_coord = np.array([[1, 0.5, 0.3], [1, 0.5, 0]])
    restri_nodes = fc.get_nodes_by_coord(coord, restri_coord)
    restri_matrix = np.empty((len(restri_nodes), 4))
    restri_matrix[:, 0] = restri_nodes
    restri_matrix[:, [1,2,3]] = np.ones((len(restri_nodes), 3), dtype='int')


The resulting mesh is shown in figure 13.

.. figure:: /conteudo/images/3d/ex4.PNG
   :scale: 50 %
   :align: center

   Figure 13: Mesh with X, Y and Z axes constrained of the nodes at coordinates: A = (1, 0.5, 0.3) and B = (1, 0.5, 0).


Constrain Z displacement of the sixth node
------------------------------

Create a matrix with the constrained node. In this example, the second and third columns are equal to 0 and the fourth is equal to 1.

.. code-block:: python

    restri_matrix = np.empty((1, 4))
    restri_matrix[0] = 6
    restri_matrix[0, [1,2]] = np.zeros((1, 2), dtype='int')
    restri_matrix[0, 3] = 1


The resulting mesh is shown in figure 14.

.. figure:: /conteudo/images/3d/ex5.PNG
   :scale: 50 %
   :align: center

   Figure 14: Mesh with Z-axis constrained of the sixth node.

