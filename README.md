# TOPDyn: Python-based Topology Optimization procedure for time-harmonic vibration problems

TOPDyn is a numerical procedure written in Python for academic purposes. The procedure is under development, focused on 2D-problems, and is based on the works:

- A critical analysis of using the dynamic compliance as objective function in topology optimization of one-material structures considering steady-state forced vibration problems.
https://www.sciencedirect.com/science/article/abs/pii/S0022460X18308563

- On the use of active and reactive input power in topology optimization of one-material structures considering steady-state forced vibration problems.
https://www.sciencedirect.com/science/article/abs/pii/S0022460X19305516

- A strategy based on the strain‐to‐kinetic energy ratio to ensure stability and convergence in topology optimization of globally resonating one‐material structures. https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6374

Some features to be implemented:

- Projection filters;

- New convergence criteria;

- Passive elements.

## Prerequisites

You will need Python 3.6.8 and the pip package manager.

It's recommended to use a [virtual environment](https://towardsdatascience.com/why-you-should-use-a-virtual-environment-for-every-python-project-c17dab3b0fd0). To do this, run the following commands:

- Linux/macos:

	Install 

		pip install virtualenv

	Create 

		virtualenv ../venv -p python3

	Activate

		source ../venv/bin/activate

- Windows:

	Install 
	
		python -m pip install virtualenv

	Create 

		virtualenv ../venv -p python3

	Activate

		../venv/Scripts/activate

To install the requirements execute

	pip install -r requirements.txt

When finished you can execute the command `deactivate` to leave the virtual environment.

Update your video card drivers.

## APIs

- SolverFEM2D: Finite element method in 2D. 

- SolverFEM3D: Finite element method in 3D (not integrated into Topology Optimization).

- Optimization algoritm: Based on MMA and GCMMA implemented by [Arjen Deetman](https://github.com/arjendeetman/GCMMA-MMA-Python). 

- More about MMA and GCMMA: [Krister Svanberg](https://people.kth.se/~krille/).

Each API has a main file with an usage example.

## Documentation 

Can be found [here](https://topdyn.readthedocs.io/en/latest/).

Need help? Have any question, or wish to ask for a missing feature? Open a issue (or send an [email](ana.rocha@mopt.com.br)).

## Authors

The authors are members of MOPT - Multidisciplinary Optimization Group, from Federal University of Santa Catarina (Florianópolis, SC, Brazil):

[Ana P. Rocha](https://www.linkedin.com/in/ana-paula-da-rocha/) and
Olavo M. Silva.

*The authors are thankful to Jacson G. Vargas and Arjen Deetman.

##
    

![alt text](https://open-pulse.github.io/OpenPulse/doc/MOPT.JPG?raw=true)

**
