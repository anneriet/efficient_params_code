
import sys
import numpy as np
import os
import scipy.optimize as opt
from math import pi
import time
import random
import networkx as nx

from pyquil import Program, get_qc
from pyquil.gates import *

# max-cut pyquil
def max_cut_gen_program(problem_size, edges, steps, qsim):
    # QAOA max-cut circuit
    prog = Program()

    params = prog.declare('angles', memory_type='REAL', memory_size=2*steps)
    ro = prog.declare('ro', memory_type='BIT', memory_size=problem_size)

    for qubit in range(problem_size):
        prog += H(qubit)
    for it in range(steps):
        for e in edges:
            prog += CNOT(e[0], e[1])
            prog += RZ(params[it], e[1])
            prog += CNOT(e[0], e[1])
        for qb in range(problem_size):
            prog += RX(params[steps+it], qb)
    for q in range(problem_size):
        prog += MEASURE(q, ro[q])
    # prog.wrap_in_numshots_loop(shots=rep)
    exe = qsim.compile(prog)
    # total_compile_time += time.time()
    return exe

# (params, graph, p, rep, qvm)
def max_cut_pyquil(program, anglelst, qvm):
    program.write_memory(region_name='angles', value=anglelst)
    open(os.path.join(curdir + "/test_output/pyquil_output.qasm"), 'w').write(str(program))

total_number_of_executions = 100

resultlst = []
for _ in range(25):
    steps = 3
    problem_size = 15
    for _ in range (total_number_of_executions):
        start_time = time.time()
        iterations = random.choice([1,2,4,6,8,10,12,14,16])

        angles = list(list(random.random() for _ in range(2*steps)) for _ in range(iterations))


        curdir = os.path.dirname(__file__)


        init_param = [0, 0] * steps
        edges = nx.random_regular_graph(2, problem_size).edges()

        qvm = get_qc('15q-qvm')
        program = max_cut_gen_program(problem_size, edges, steps, qvm)
        for i in range(iterations):
            max_cut_pyquil(program, angles[i], qvm)
         total_time = time.time()  - start_time # end_time - start_time - (end_sim - start_sim)

        total_compile_time = total_time
        total_simulation_time = 0
        filename = 'pyquil_wo_optimiser.txt'
        curdir = os.path.dirname(__file__)

        line = (time.strftime("%d/%m/%y %H:%M:%S") + "\t" +  # +"Total time:\t"
                                                str(total_time) + "\t" +  #Compile time: \t" + 
                                                str(total_compile_time) +"\t" + # Simulation time:\t" + 
                                                str(total_simulation_time) + "\t" +  # problem_size:\t" + 
                                                str(problem_size) + "\t" +  #steps: \t" + 
                                                str(steps) + "\t" + # iterations:\t" + 
                                                str(iterations) + '\n')
        # open a file to append
        with open(os.path.join(curdir, filename), 'a') as f:
            f.write(line)
        print(line)
