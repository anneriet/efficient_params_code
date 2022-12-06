import numpy as np
from qiskit import QuantumCircuit, transpile, execute, assemble
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
from qiskit.circuit import Parameter, ParameterVector


import time
import networkx as nx

import scipy.optimize as opt
import random as rn
import os
# max-cut
def eval_cost_qiskit(resultlst, edges):
    # evaluate Max-Cut
    c = 0
    for res, count in resultlst:
        for e in edges:
            c += count * (0.5 * (1 - int(res[e[0]]) * int(res[e[1]])))
    return c 


# return a regular graph
def regular_graph(n):
    edges = []
    for i in range(n-1):
        edges.append([i, i+1])
    edges.append([0, n-1])
    for i in range(n-2):
        edges.append([i, i+2])
    edges.append([0, n-2])
    edges.append([1, n-1])
    return edges

def max_cut_circuit(problem_size, edges, steps, init_param):
    global total_compile_time
    total_compile_time -= time.time()

    qc = QuantumCircuit(problem_size, problem_size)

    gammas = ParameterVector("gammas", steps)
    betas = ParameterVector("betas", steps)
    
    for q in range(problem_size):
        qc.h(q)
        
    for it in range(steps):
        for e in edges:
            qc.cnot(e[0], e[1])
            qc.rz(gammas[it], e[1])
            qc.cnot(e[0], e[1])
        for qb in range(problem_size):
            qc.rx(betas[it], qb)
    for q in range(problem_size):
        qc.measure(q, q)
    total_compile_time += time.time()
    return [qc, gammas, betas]

def max_cut_qiskit(qc, anglelst, edges, steps, rep, qsim):
    global total_compile_time
    global total_simulation_time
    beta = anglelst[0:steps]
    gamma = anglelst[steps:2*steps]
    for pb in range(steps):
        if beta[pb] < 0 or beta[pb] > (2 * np.pi) or gamma[pb] < 0 or gamma[pb] > (2 * np.pi):
            return 0

    total_compile_time -= time.time()
    compiled_circuit = assemble([qc.bind_parameters({gammas : gamma, betas : beta})], qsim)
    total_compile_time += time.time()
    total_simulation_time -= time.time()
    # Execute the circuit on the qasm simulator
    job = qsim.run(compiled_circuit, shots=rep)
    # Grab results from the job
    result = job.result()
    # Returns counts
    counts = result.get_counts().items()
    # print("data: ", counts)
    total_simulation_time += time.time()
    return eval_cost_qiskit(counts, edges)


def mcp_inv_params(angles, transpiled_circuit, edges, steps,rep, qsim):
    out = max_cut_qiskit(transpiled_circuit, angles, edges, steps, rep, qsim)
    if out == 0:
        inv = 2  
    else:
        inv = 1 / out
    return inv


def nm_params(init_param, transpiled_circuit, edges, steps, rep, qsim, maxevals):
    
    res = opt.minimize(mcp_inv_params, init_param, args=(transpiled_circuit, edges, steps, rep, qsim), method='nelder-mead',
                       options={'fatol': 1e-5, 'maxfev': maxevals, 'disp': False})
    return [1 / res.fun, res.x, res.success, res.nit, res.nfev]


problem_size_lst = [3,4,5,6,7,8] 
for _ in range(100):
    steps = 3
    problem_size = rn.choice(problem_size_lst)
    opt_lvl = 1
    total_simulation_time = 0
    total_compile_time = 0
    start_time = time.time()


    rep = 10
    maxevals = 100



    init_param = [0, 0] * steps
    edges = nx.random_regular_graph(2, problem_size).edges()
    simulator = QasmSimulator()

    qiskit_circuit, gammas, betas = max_cut_circuit(problem_size, edges, steps, init_param)
    transpiled_circuit = transpile(qiskit_circuit, simulator, optimization_level=opt_lvl)
    result = nm_params(init_param, transpiled_circuit, edges, steps, rep, simulator, maxevals)
                
    total_time = time.time()  - start_time # end_time - start_time - (end_sim - start_sim)

    filename = 'qiskit_maxcut.txt'
    curdir = os.path.dirname(__file__)

    line = (time.strftime("%d/%m/%y %H:%M:%S") + "\t" +  # +"Total time:\t"
                                            str(total_time) + "\t" +  #Compile time: \t" + 
                                            str(total_compile_time) +"\t" + # Simulation time:\t" + 
                                            str(total_simulation_time) + "\t" +  # problem_size:\t" + 
                                            str(problem_size) + "\t" +  #steps: \t" + 
                                            str(steps) + "\t" + # maxiter:\t" + 
                                            str(maxevals) + "\t" +  #nit:\t" + 
                                            str(result[3]) + "\t" +  #nfev:\t" + 
                                            str(result[4]) + "\t" + # opt_level:\t
                                            str(opt_lvl) + "\n" )
    # open a file to append
    with open(os.path.join(curdir, filename), 'a') as f:
        f.write(line)
