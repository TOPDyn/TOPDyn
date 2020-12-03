Constrain nodes displacements in 2D
=======================================

Figure 6 shows a mesh with *nelx = 4*, *nely = 3*, and how the nodes are numbered. The nodes are numbered from left to right and from bottom to top.

.. figure:: /conteudo/images/2d/fig7.png
   :scale: 50 %
   :align: center

   Figure 6: Mesh 2D with numbered nodes.


Constrain X displacement of the nodes with Y = 0.2
-------------------------------------------------------

It is necessary to find the nodes with Y = 0.2 with function :meth:`SolverFEM2D.functions2D.get_nodes1d`. Then a matrix with three columns is created with the constrained nodes. The first column is the node number, the second column designates whether the node is constrained on X displacement (*1: True*, *0: False*), and the third column designates whether the node is constrained on Y displacement (*1: True*, *0: False*). In this example, the second column is equal to 1, and the third is equal to 0.

.. code-block:: python

    restri_nodes = fc.get_nodes1d(coord, 0.2, 0.001, 2)
    restri_matrix = np.empty((len(restri_nodes), 3), dtype='int')
    restri_matrix[:, 0] = restri_nodes
    restri_matrix[:, 1] = np.ones(restri_nodes.shape[0])
    restri_matrix[:, 2] = np.zeros(restri_nodes.shape[0])

The resulting mesh is shown in figure 7.

.. figure:: /conteudo/images/2d/fig1.png
   :scale: 50 %
   :align: center

   Figure 7: Mesh with X-axis nodes constrained with Y = 0.2.

Constrain X and Y displacements of the nodes with X = 0
-------------------------------------------------------------

The objective is to obtain the nodes whose X-axis (column = 1) is equal to 0 with a margin :math:`\epsilon` = 0.001. The second and third columns of the constrained nodes matrix are equal to 1.

.. code-block:: python

    restri_nodes = fc.get_nodes1d(coord, 0, 0.001, 1)
    restri_matrix = np.empty((len(restri_nodes), 3), dtype='int')
    restri_matrix[:, 0] = restri_nodes
    restri_matrix[:, [1, 2]] = np.ones((restri_nodes.shape[0], 2))

The resulting mesh is shown in figure 8.

.. figure:: /conteudo/images/2d/fig3.png
   :scale: 50 %
   :align: center

   Figure 8: Mesh with all X and Y axes elements constrained with X = 0.

Constrain X and Y displacements of the node at coordinate (1, 0.23)
----------------------------------------------------------------------

This example shows how the rounding of function :meth:`SolverFEM2D.functions2D.get_nodes_by_coord` works. Unique values on the Y-axis don't include 0.23:

.. code-block:: python

    array([0.  , 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 ])

The algorithm will round to 0.25, and the node obtained is at coordinate (1, 0.25). Besides, the second and third columns of the constrained nodes matrix are equal to 1.

.. code-block:: python

    array([[66.  ,  1.  ,  0.25]])

.. code-block:: python

    restri_coord = np.array([[1, 0.23]])
    restri_nodes = fc.get_nodes_by_coord(coord, restri_coord)
    restri_matrix = np.empty((len(restri_nodes), 3), dtype='int')
    restri_matrix[:, 0] = restri_nodes
    restri_matrix[:, [1, 2]] = np.ones((restri_nodes.shape[0], 2))

The resulting mesh is shown in figure 9.

.. figure:: /conteudo/images/2d/fig4.png
   :scale: 50 %
   :align: center

   Figure 9: Mesh with X and Y axes constrained of the node at coordinate (1, 0.23).

Constrain the eighth node displacement
------------------------------------------

Create a matrix with the constrained node. In this example, the second and third columns are equal to 1.

.. code-block:: python

    restri_matrix = np.empty((1, 3), dtype='int')
    restri_matrix[:, 0] = 8
    restri_matrix[:, [1, 2]] = np.ones((1, 2))


The resulting mesh is shown in figure 10.


.. figure:: /conteudo/images/2d/fig6.png
   :scale: 50 %
   :align: center

   Figure 10: Mesh with X and Y axes constrained of the eighth node.





