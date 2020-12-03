These how-to guides will step you through common tasks in using and configuring a TOPDyn setup.

Figure 1 shows a diagram with the structure of TOPDyn.

.. figure:: /conteudo/images/diagram.png
   :scale: 50 %
   :align: center

   Figure 1: TOPDyn structure.

Box constraint to Area and R ratio
================
An interval is defined for the Area value or R Ratio. For that, put a positive value that is the upper limit and a negative value that is the lower limit. Figure 2 reproduces the constrain of the code below.

.. code-block:: python

   constr_func = ['Area', 'Area']
   constr_values = [30, -10]

.. figure:: /conteudo/images/box_constrain.png
   :scale: 50 %
   :align: center

   Figure 2: Area constraint.


