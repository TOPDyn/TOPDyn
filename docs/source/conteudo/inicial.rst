Constraint functions
======================================

It's possible to define a maximum value for the constraint or a box constraint like showed in the Figure 2.

The constraint functions that can be used are:
   - Area
   - Strain-to-kinetic energy ratio
   - Compliance
   - Local kinetic energy
   - Local elastic potential energy
   - Local strain-to-kinetic energy ratio

In the code they are passed respectively as:
   - Area -> area
   - Strain-to-kinetic energy ratio -> r_ratio
   - Compliance -> compliance
   - Local kinetic energy -> local_ki
   - Local elastic potential energy -> local_ep
   - Local strain-to-kinetic energy ratio -> local_r

For the last four constraint functions (compliance, local_ki, local_ep, local_r) is necessary to pass the frequency that will be used. For example:

.. code-block:: python

   constr_func = ['local_ep']
   constr_values = [(70, 50)]

It means that the constraint for the local elastic potential energy is 70 and the frequency that it will be calculate is 50 Hz.

To use a box constraint is necessary to put a positive value that is the upper limit and a negative value that is the lower limit. Figure 1 reproduces the constrain of the code below.

.. code-block:: python

   constr_func = ['Area', 'Area']
   constr_values = [30, -10]

.. figure:: /conteudo/images/box_constrain.png
   :scale: 70 %
   :align: center

   Figure 1: Area constraint.


