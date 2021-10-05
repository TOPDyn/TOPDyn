Add loads in 2D
===================


It's passed a list of dictionaries. There are two options for add loads:

   - Add loads by column 

   - Add loads by coordinate


.. Hint::

   It's possible combine and use both options.


Add loads by coordinate
--------------------------

The dictionary has the following keys: 

   - "x_coord": coordinate in X-axis.
   
   - "y_coord": coordinate in Y-axis.
   
   - "x_direc": indicates if is applied a load in X-direction.
         - If 1 indicates that the load is applied in the negative direction of the X-direction.
         - If 1 indicates that the load is applied in the positive direction of the X-direction.
         - If 0 indicates that no load is applied in the X-direction.
   
   - "y_direc": indicates if is applied a load in Y-direction.
         - If 1 indicates that the load is applied in the negative direction of the Y-direction.
         - If 1 indicates that the load is applied in the positive direction of the Y-direction.
         - If 0 indicates that no load is applied in the Y-direction.
   
   - "force": The module of the load in Newton.

In the Figure 3 two loads are applied in the mesh:

   - A negative force of modulus 100 N is applied in the Y-direction to the node at coordinate (1,0.5).

   - A positive force of modulus 100 N is applied in the X-direction to the node at coordinate (0,0).

In Python language:

.. code-block:: python

   load_matrix = [{"x_coord":1, "y_coord":0.5, "x_direc":0, "y_direc":-1, "force":100}, \
                   {"x_coord":0, "y_coord":0, "x_direc":1, "y_direc":0, "force":100}]


.. figure:: /conteudo/images/force_node.png
   :scale: 50 %
   :align: center

   Figure 3: 2D mesh with load applied at first and last node.


Add loads by column
-------------------------

The dictionary has the following keys: 

   - "coord": Coordinate.
   
   - "axis": Direction that the load is applied. It can be: 1 (X-Axis) or 2 (Y-Axis).

   - "eps": Margin of error.
   
   - "x_direc": indicates if is applied a load in X-direction.
         - If 1 indicates that the load is applied in the negative direction of the X-direction.
         - If 1 indicates that the load is applied in the positive direction of the X-direction.
         - If 0 indicates that no load is applied in the X-direction.
   
   - "y_direc": indicates if is applied a load in Y-direction.
         - If 1 indicates that the load is applied in the negative direction of the Y-direction.
         - If 1 indicates that the load is applied in the positive direction of the Y-direction.
         - If 0 indicates that no load is applied in the Y-direction.
   
   - "force": The module of the load in Newton.

In the example it's applied a positive force of modulus 100 N in X-direction to all the nodes with Y=0.1.

.. code-block:: python

   load_matrix = [{"coord":0.1, "axis":2, "x_direc":1, "y_direc":0, "force":100, "eps":0.001}]

The Figure 4 shows the load applied in the mesh as the example.

.. figure:: /conteudo/images/force_all.png
   :scale: 70 %
   :align: center

   Figure 4: Load applied at all nodes with Y = 0.1 in 2D mesh.