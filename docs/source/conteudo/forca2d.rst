Add loads in 2D
================


Add loads by node number
----------------

Create a matrix with the loads. The first column is the node number, second and third columns indicate how the load is applied in X and Y axes, respectively. The values can be:

   - -1 indicates that the load is applied in the negative direction of the axis.
    
   - 1 indicates that the load is applied in the positive direction of the axis.
    
   - 0 indicates that no load is applied in the axis.

The last column (fourth) is the module of the load in Newton. Each row of the matrix is a different load. In this example, the first node was loaded at the X-axis in the positive direction, and the last node has loaded at Y-axis in the negative direction. 

.. code-block:: python

   force_matrix = np.empty((2, 4))
   force_matrix[:, 0] = [1, (nelx + 1) * (nely + 1)]
   force_matrix[:, [1,2,3]] = [[1, 0, 10], [0, -1, 100]]

The resulting mesh is shown in figure 3.

.. figure:: /conteudo/images/force_node.png
   :scale: 50 %
   :align: center

   Figure 3: Mesh with load applied at firts and last node.


Add loads by coordinate
------------------

Firstly, it's necessary to find the nodes which correspond to the coordinates provided. Function :meth:`SolverFEM2D.functions2D.get_nodes_by_coord` obtains the coordinate nodes passed in the second argument that are the same as the coordinate matrix generate by function :meth:`SolverFEM2D.functions2D.regularmeshQ4`. 

Each line of array *restri_coord* is a coordinate to compare with the coordinate matrix, and the columns are the X and Y coordinates, respectively. If the provided coordinate doesn't exist in the coordinate matrix. Then the algorithm will round to the nearest.

In this example, will be added a load in Y-axis in the negative direction at the coordinate A = (1, 0) and a load in X-axis in the positive direction B = (1, 0.5). 

.. code-block:: python

   restri_coord = np.array([[1, 0], [1, 0.5]])
   nodes_apply_F = fc.get_nodes_by_coord(coord, restri_coord)
   force_matrix = np.empty((2, 4))
   force_matrix[:, 0] = nodes_apply_F
   force_matrix[:, [1,2,3]] = [[0, -1, 10], [1, 0, 100]]

The resulting mesh is shown in figure 4.

.. figure:: /conteudo/images/force_coord.png
   :scale: 50 %
   :align: center

   Figure 4: Mesh with load applied by coordinate.


Add loads in nodes with Y = 0.1
------------------

Function :meth:`SolverFEM2D.functions2D.get_nodes1d` obtains the nodes of the coordinate matrix whose Y-axis (column = 2) is equal to 0.1 with a margin :math:`\epsilon` (in the code is 0.001). The more discretized the mesh, the smaller the acceptable :math:`\epsilon` margin must be.

After creating the load matrix, the objective is to apply force in the positive X direction, so the second column (column 1) will be filled with the value 1.


.. code-block:: python

    nodes_apply_F = fc.get_nodes1d(coord, 0.1, 0.001, 2)
    force_matrix = np.empty((nodes_apply_F.shape[0], 4))
    force_matrix[:, 0] = nodes_apply_F
    force_matrix[:, 1] = np.ones(nodes_apply_F.shape[0])
    force_matrix[:, 2] = np.zeros(nodes_apply_F.shape[0])
    force_matrix[:, 3] = 10 * np.ones(nodes_apply_F.shape[0])

The resulting mesh is shown in figure 5.

.. figure:: /conteudo/images/force_all.png
   :scale: 50 %
   :align: center

   Figure 5: Load applied in all nodes with Y = 0.1.


