Constrain nodes displacements in 3D
=======================================

Figure 10 shows a mesh with *nelx = 2*, *nely = 2*, *nelz = 1*, and how the nodes are numbered. 

.. figure:: /conteudo/images/3d/mesh3d.JPG
   :scale: 65 %
   :align: center

   Figure 10: 3D mesh with numbered nodes.

It's passed a list of dictionaries. There are two options for constrain nodes displacements:

   - Constrain nodes displacements by column 

   - Constrain nodes displacements by coordinate

.. Hint::

   It's possible combine and use both options.


Constrain nodes displacements by coordinate
----------------------------------------------------

The dictionary has the following keys: 

   - "x_coord": coordinate in X-axis.
   
   - "y_coord": coordinate in Y-axis.

   - "z_coord": coordinate in Z-axis.
   
   - "constrain_disp_x": indicates if the node displacement is contraint in X-direction.
         - If 1 constrain the node displacement in the X-direction.
         - If 0 don't constrain the node displacement in the X-direction.

    - "constrain_disp_y": indicates if the node displacement is contraint in Y-direction.
         - If 1 constrain the node displacement in the Y-direction.
         - If 0 don't constrain the node displacement in the Y-direction.

    - "constrain_disp_z": indicates if the node displacement is contraint in Z-direction.
         - If 1 constrain the node displacement in the Z-direction.
         - If 0 don't constrain the node displacement in the Z-direction.

In the example is constrained the node displacement at coordinateS (1, 0.5, 0.3) and (1, 0.5, 0).

.. code-block:: python
    
    restri_matrix = [{"x_coord":1, "y_coord":0.5, "z_coord":0.3, "constrain_disp_x":1, "constrain_disp_y":1, "constrain_disp_z":1}, \
                     {"x_coord":1, "y_coord":0.5, "z_coord":0, "constrain_disp_x":1, "constrain_disp_y":1, "constrain_disp_z":1}]

The resulting mesh is shown in figure 11.

.. figure:: /conteudo/images/3d/ex4.PNG
   :scale: 50 %
   :align: center

   Figure 11: Mesh with nodes at coordinates (1, 0.5, 0.3) and (1, 0.5, 0).


Constrain nodes displacements by column
----------------------------------------------------

The dictionary has the following keys: 

   - "coord": coordinate.
   
   - "axis": Direction that the node displacement is constraint. It can be: 1 (X-Axis), 2 (Y-Axis) or 3 (Z-Axis).

   - "eps": Margin of error.
   
   - "constrain_disp_x": indicates if the node displacement is contraint in X-direction.
         - If 1 constrain the node displacement in the X-direction.
         - If 0 don't constrain the node displacement in the X-direction.

    - "constrain_disp_y": indicates if the node displacement is contraint in Y-direction.
         - If 1 constrain the node displacement in the Y-direction.
         - If 0 don't constrain the node displacement in the Y-direction.

    - "constrain_disp_z": indicates if the node displacement is contraint in Z-direction.
         - If 1 constrain the node displacement in the Z-direction.
         - If 0 don't constrain the node displacement in the Z-direction.

In the example is constrained the X displacement of the nodes in the plane Z = 0.3.

.. code-block:: python

    restri_matrix = [{"coord":0.3, "axis":3, "eps":0.001, "constrain_disp_x":1, "constrain_disp_y":0, "constrain_disp_z":0}]

The resulting mesh is shown in figure 12.

.. figure:: /conteudo/images/3d/ex2.PNG
   :scale: 50 %
   :align: center

   Figure 12: Nodes in the plane Z = 0.3.