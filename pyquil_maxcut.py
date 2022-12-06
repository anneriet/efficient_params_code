
import sys

import numpy as np
import os
import scipy.optimize as opt
from math import pi
import time
import random as rn
import networkx as nx


from pyquil import Program, get_qc
from pyquil.gates import *

def eval_cost_pyquil(vallst, edges):
    # evaluate Max-Cut
    global resultlst
    c = 0
    for res in vallst:
        for e in edges:
            c += 0.5 * (1 - int(res[e[0]]) * int(res[e[1]]))
    resultlst.append(c)
    return c

def max_cut_gen_program(problem_size, edges, steps, init_param, rep, qsim):
    global total_compile_time
    total_compile_time -= time.time()
    # QAOA max-cut circuit
    prog = Program()

    gammas = prog.declare('gamma', memory_type='REAL', memory_size=steps)
    betas = prog.declare('beta', memory_type='REAL', memory_size=steps)
    ro = prog.declare('ro', memory_type='BIT', memory_size=problem_size)

    for qubit in range(problem_size):
        prog += H(qubit)
    for it in range(steps):
        for e in edges:
            prog += CNOT(e[0], e[1])
            prog += RZ(-gammas[it], e[1])
            prog += CNOT(e[0], e[1])
        for qb in range(problem_size):
            prog += RX(2*betas[it], qb)
    for q in range(problem_size):
        prog += MEASURE(q, ro[q])
    prog.wrap_in_numshots_loop(shots=rep)
    exe = qsim.compile(prog)
    total_compile_time += time.time()
    return exe

def max_cut_pyquil(program, anglelst, edges, steps, rep, qsim):
    global total_compile_time
    global total_simulation_time
    beta = anglelst[0:steps]
    gamma = anglelst[steps:2*steps]
    for pb in range(steps):
        if beta[pb] < 0 or beta[pb] > (2 * np.pi) or gamma[pb] < 0 or gamma[pb] > (2 * np.pi):
            return 0

    program.write_memory(region_name='gamma', value=gamma)
    program.write_memory(region_name='beta', value=beta)
    total_simulation_time -= time.time()
    results = qsim.run(program)
    total_simulation_time += time.time()
    return eval_cost_pyquil(results.readout_data.get('ro'), edges)


def mcp_inv_pyquil(angles, program, edges, steps,rep, qsim):
    out = max_cut_pyquil(program, angles, edges, steps, rep, qsim)
    if out == 0:
        inv = 2 
    else:
        inv = 1 / out
    return inv

def nm_pyquil(q_func_str, init_param, program, edges, steps, rep, qsim, maxevals):
    min_func = {
        'mcp': mcp_inv_pyquil
    }.get(q_func)

    res = opt.minimize(min_func, init_param, args=(program, edges, steps, rep, qsim), method='nelder-mead',
                       options={'fatol': 1e-5, 'maxfev': maxevals, 'disp': False})
    return [1 / res.fun, res.x, res.success, res.nit, res.nfev]

problem_size_lst = [3,4,5,6,7,8] # program length

total_number_of_executions = 1

resultlst = []
for _ in range(25):
    steps = 3
    problem_size = rn.choice(problem_size_lst)
    for _ in range (total_number_of_executions):
        total_simulation_time = 0
        total_compile_time = 0
        start_time = time.time()


        rep = 10
        maxevals = 100



        init_param = [0, 0] * steps
        edges = nx.random_regular_graph(2, problem_size).edges()
        q_func = 'mcp'
        qvm = get_qc('10q-qvm')
        program = max_cut_gen_program(problem_size, edges, steps, init_param, rep, qvm)
        result = nm_pyquil(q_func, init_param, program, edges, steps, rep, qvm, maxevals)
        total_time = time.time()  - start_time # end_time - start_time - (end_sim - start_sim)

        filename = 'pyquil_maxcut.txt'
        curdir = os.path.dirname(__file__)

        line = (time.strftime("%d/%m/%y %H:%M:%S") + "\t" +  # +"Total time:\t"
                                                str(total_time) + "\t" +  #Compile time: \t" + 
                                                str(total_compile_time) +"\t" + # Simulation time:\t" + 
                                                str(total_simulation_time) + "\t" +  # problem_size:\t" + 
                                                str(problem_size) + "\t" +  #steps: \t" + 
                                                str(steps) + "\t" + # maxiter:\t" + 
                                                str(maxevals) + "\t" +  #nit:\t" + 
                                                str(result[3]) + "\t" +  #nfev:\t" + 
                                                str(result[4]) + "\n" )
        # open a file to append
        with open(os.path.join(curdir, filename), 'a') as f:
            f.write(line)


# plt.plot(resultlst)
# plt.show()