# Efficient parameterised compilation for hybrid quantum programming
This repository contains the code that was used to generate the results for the paper https://arxiv.org/abs/2208.07683. In order to run the experiments, please first install the application discussed in the Setup section below.

### Setup
[OpenQL_PC and OpenQL](https://github.com/anneriet/OpenQL/tree/parameters_for_hybrid_comp), follow the instructions for installing from source in the [OpenQL documentation](https://openql.readthedocs.io/).    
Qiskit (0.37.0): ``pip install qiskit==0.3.7``  
PyQuil (3.1.0): ``pip install pyquil==3.1.0``  
And install the qvm and qvc from [here](https://pyquil-docs.rigetti.com/en/stable/start.html#upgrading-or-installing-pyquil).  

### Timing with random angles
These files were used for the timing of the MAXCUT circuit but with randomised angles as parameters for each iteration, so that the measured compile time does not include the classical optimiser (which would otherwise generate the angles).  
params_wo_optimiser.py: Timing for OpenQL, OpenQL_PC and Qiskit.  
pyquil_wo_optimiser.py: Timing for PyQuil. PyQuil compilation requires running the qvc compiler as a separate process. To get accurate results the PyQuil tests were not run concurrently with OpenQL, OpenQL_PC and Qiskit.  

### Full MAXCUT
Code used for the timing of the full MAXCUT algorithm, including simulation and the classical optimiser:  
openql_openql_pc_maxcut.py  
qiskit_maxcut.py  
pyquil_maxcut.py  
Separate files were used because of conflicting code dependencies in the simulators.
