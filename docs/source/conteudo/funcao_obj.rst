Objective functions 
=====================================


TOPDyn has the following objective functions: 

- Input power

- Compliance

- Elastic Potential Energy

- Kinetic Energy

- R Ratio


To use one objective function is necessary to pass to the argument **func_name** the name of the function as a string, the frequency that pretends to be optimized, argument **freq1** (the default value is 180 Hz), and make sure the weight applied (parameter **n1**) to function is 1. If the value of n1 is negative, then the function is maximized. The code below shows how is the definition of these parameters in Python language.

.. code-block:: python

    func_name = 'Elastic Potential Energy'
    freq1 = 275
    n1 = 1


In order to be able to optimize two objective functions, they must be combined to minimize both. Thus, the following multi-objective function is proposed:


.. math::

   F = n1 \cdot F1 + (1 - abs(n1)) \cdot F2

where F1 is the function passed in the parameter **func_name**, and F2 is the function passed in the parameter **multiobjective**, and **n1** is the weight that controls the priority between the functions. The parameter **multiobjetive** is a tuple with the name of the function and the frequency associated with this function. In the code example below, the value of n1 is 0.99, so the input power function will have a higher priority than the compliance function.

.. code-block:: python

    multiobjective = ('Compliance', 0)
    func_name = 'Input Power'
    freq1 = 245
    n1 = 0.99

