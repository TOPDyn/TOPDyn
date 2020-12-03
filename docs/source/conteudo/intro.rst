Features
================

- SolverFEM2D: Finite element method in 2D. 

- SolverFEM3D: Finite element method in 3D. 

- Optimization algoritm: Based on MMA and GCMMA implemented by `Arjen Deetman <https://github.com/arjendeetman/GCMMA-MMA-Python>`_. 

- More about MMA and GCMMA: `Krister Svanberg <https://people.kth.se/~krille/>`_.

Each feature has a main file with an usage example.


Prerequisites
================

You will need python 3.6.8 and the pip package manager.

It's recommended to use a `virtual environment <https://towardsdatascience.com/why-you-should-use-a-virtual-environment-for-every-python-project-c17dab3b0fd0>`_. To do this, run the following commands:

- Linux/macos:

	Install::

		$ pip install virtualenv

	Create::

		$ virtualenv ../venv -p python3

	Activate::

		$ source ../venv/bin/activate

- Windows:

	Install:: 
	
		$ python -m pip install virtualenv

	Create:: 

		$ virtualenv ../venv -p python3

	Activate::

		$ ../venv/Scripts/activate

To install the requirements execute:

.. code:: python

	pip install -r requirements.txt

When finished you can execute the command `deactivate` to leave the virtual environment.

Update your video card drivers.

