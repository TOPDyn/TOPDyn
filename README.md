# TOPDyn: Python-based Topology Optimization procedure for time-harmonic vibration problems

TOPDyn is a numerical procedure written in Python for academic purposes. The procedure is under development, focused on 2D-problems, and is based on the works:

- A critical analysis of using the dynamic compliance as objective function in topology optimization of one-material structures considering steady-state forced vibration problems.
https://www.sciencedirect.com/science/article/abs/pii/S0022460X18308563

- On the use of active and reactive input power in topology optimization of one-material structures considering steady-state forced vibration problems.
https://www.sciencedirect.com/science/article/abs/pii/S0022460X19305516

- A strategy based on the strain‐to‐kinetic energy ratio to ensure stability and convergence in topology optimization of globally resonating one‐material structures. https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6374

- Shape and position preserving design of vibrating structures by controlling local energies through 
topology optimization. https://www.sciencedirect.com/science/article/pii/S0022460X21005150

Feature to be implemented:

- Projection filters;

## Prerequisites

You will need Python 3.6.8 and the pip package manager.

To install the requirements execute

	pip install -r requirements.txt

Update your video card drivers.

## APIs

- SolverFEM2D: Finite element method in 2D. 

- SolverFEM3D: Finite element method in 3D (not integrated into Topology Optimization).

- Optimization algoritm: Based on MMA and GCMMA implemented by [Arjen Deetman](https://github.com/arjendeetman/GCMMA-MMA-Python). 

- More about MMA and GCMMA: [Krister Svanberg](https://people.kth.se/~krille/).

> :raised_hand: Each API has a **main file** with an usage example and It's the place where you must put your code. :raised_hand:

## :closed_book: Documentation 

Can be found [here](https://topdyn.readthedocs.io/en/latest/).

Need help? Have any question, or wish to ask for a missing feature? Open a issue.

## Authors

The authors are members of MOPT - Multidisciplinary Optimization Group, from Federal University of Santa Catarina (Florianópolis, SC, Brazil):

[Ana P. Rocha](https://www.linkedin.com/in/ana-paula-da-rocha/) and
Olavo M. Silva.

*The authors are thankful to Jacson G. Vargas and Arjen Deetman.

##
    

![alt text](https://open-pulse.github.io/OpenPulse/doc/MOPT.JPG?raw=true)

**
