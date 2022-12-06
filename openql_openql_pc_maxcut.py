#!/usr/bin/python3

import sys

from openql import openql as ql
from qxelarator import qxelarator
import numpy as np
import os
import scipy.optimize as opt
from math import pi
import time
import random as rn
import networkx as nx

# return a regular graph
def regular_graph(n):
    edges = []
    for i in range(n-1):
        edges.append([i, i+1])
    edges.append([0, n-1])
    if n > 4:
        for i in range(n-2):
            edges.append([i, i+2])
        edges.append([0, n-2])
    edges.append([1, n-1])
    return edges

# max-cut
def eval_cost_openql(res, edges):
    # evaluate Max-Cut
    c = 0
    for e in edges:
        c += 0.5 * (1 -res[e[0]] * res[e[1]])
    return c 

def max_cut_gen_kernel(k, problem_size, edges, steps, init_param):
    global total_compile_time
    total_compile_time -= time.time()
    ql_paramlst = []

    for init_beta in init_param[0:steps]:
        ql_paramlst.append(ql.Param("ANGLE", 2*init_beta)) # beta

    for init_gamma in init_param[steps:2*steps]: # DOES THIS MATTER
        ql_paramlst.append(ql.Param("ANGLE", -init_gamma)) # gamma

    for q in range(problem_size):
        k.gate("hadamard", [q])
        
    for it in range(steps):
        for e in edges:
            k.gate("cnot", e)
            k.gate("rz", [e[1]], 0, ql_paramlst[steps+it]) # gamma
            k.gate("cnot",e)
        for qb in range(problem_size):
            k.gate("rx", [qb], 0, ql_paramlst[it]) # beta
    # for q in range(problem_size):
    #     k.gate("measure", [q])
    total_compile_time += time.time()
    return [k, ql_paramlst]
    

def max_cut_openql_params(program, ql_paramlst, anglelst, edges, steps, rep, qsim):
    global total_compile_time
    global total_simulation_time
    beta = anglelst[0:steps]
    gamma = anglelst[steps:2*steps]
    for pb in range(steps):
        if beta[pb] < 0 or beta[pb] > (2 * np.pi) or gamma[pb] < 0 or gamma[pb] > (2 * np.pi):
            return 0
    angles = np.zeros(2*steps)
    angles[0:steps] = beta*2
    angles[steps:2*steps] = gamma*-1
    total_compile_time -= time.time()
    c.compile(program, ql_paramlst, angles) 
    total_compile_time += time.time()
    total_simulation_time -= time.time()
    res = qsim.execute_and_get_average_measurement(rep)
    total_simulation_time += time.time()
    return eval_cost_openql(res, edges)

def mcp_inv_params(angles, program, ql_paramlst, edges, steps,rep, qsim):
    out = max_cut_openql_params(program, ql_paramlst, angles, edges, steps, rep, qsim)
    if out == 0:
        inv = 2  
    else:
        inv = 1 / out
    return inv

def nm_params(q_func_str, init_param, program, ql_paramlst, edges, steps, rep, qsim, maxevals):
    res = opt.minimize(mcp_inv_params, init_param, args=(program, ql_paramlst, edges, steps, rep, qsim), method='nelder-mead',
                       options={'ftol': 1e-5, 'maxfev': maxevals, 'disp': False}) # is fatol something different?
    return [1 / res.fun, res.x, res.success, res.nit, res.nfev]

    
def max_cut_noparams(angles, problem_size, edges, steps, rep, qsim):
    global total_simulation_time
    global total_compile_time
    beta = angles[0:steps]
    gamma = angles[steps:2*steps]
    for pb in range(steps):
        if beta[pb] < 0 or beta[pb] > (2 * np.pi) or gamma[pb] < 0 or gamma[pb] > (2 * np.pi):
            return 0
    total_compile_time -= time.time()
    program = ql.Program(programname, platf, problem_size, 0, 0)
    k = ql.Kernel ("qk", platf, problem_size) 
    for q in range(problem_size):
        k.gate("hadamard", [q])
    for it in range(steps):
        for e in edges:
            k.gate("cnot", e)
            k.gate("rz", [e[1]], 0, -gamma[it])
            k.gate("cnot",e)
        for qb in range(problem_size):
            k.gate("rx", [qb], 0, 2*beta[it])
    program.add_kernel(k)
    c.compile(program)
    total_compile_time += time.time()
    total_simulation_time -= time.time()
    res = qsim.execute_and_get_average_measurement(rep)
    total_simulation_time += time.time()
    return eval_cost_openql(res, edges)
    
def mcp_inv_noparams(angles, problem_size, edges, steps, rep, qsim):
    out = max_cut_noparams(angles, problem_size, edges, steps, rep, qsim)
    if out == 0:
        inv = 100 
    else:
        inv = 1 / out
    return inv

def nm_noparams(q_func, init_param, problem_size, edges, steps, rep, qsim, maxevals):
    res = opt.minimize(mcp_inv_noparams, init_param, args=(problem_size, edges, steps, rep, qsim), method='nelder-mead',
                       options={'ftol': 1e-5, 'maxfev': maxevals, 'disp': False})
    return [1 / res.fun, res.x, res.success, res.nit, res.nfev]

    

problem_size_lst = [3,4,5,6,7,8] 

for _ in range(199):
    steps = 3
    problem_size = rn.choice(problem_size_lst)
    withparams = rn.randint(0,1)
    total_simulation_time = 0
    total_compile_time = 0
    start_time = time.time()

    ql.set_option('log_level', 'LOG_ERROR')
    curdir = os.path.dirname(__file__)
    output_dir = os.path.join(curdir, 'test_output')
    ql.set_option('output_dir', output_dir)
    config_fn = os.path.join(curdir, 'hardware_config_qx.json')
    c = ql.Compiler("testCompiler", config_fn)
    c.add_pass_alias("Writer", "outputIR")
    platf = ql.Platform("myPlatform", config_fn)

    programname = "maxcut_execution_time"


    rep = 10
    maxevals = 100
    init_param = list(rn.random()*2*np.pi for _ in range(2*steps))
    q_func = 'mcp'
    edges = regular_graph(problem_size)
    qx = qxelarator.QX()
    qx.set(os.path.join(output_dir + "/" + programname + ".qasm"))
    print("withparams? ", withparams, "problem_size: ", problem_size)
    filename = 'maxcut_with_without_parameters_runtime_new_03.txt'
    if(withparams):
        program = ql.Program(programname, platf, problem_size, 0, 0)
        k = ql.Kernel ("qk", platf, problem_size) 
        k, paramlst = max_cut_gen_kernel(k, problem_size, edges, steps, init_param)
        program.add_kernel(k)
        result = nm_params(q_func, init_param, program, paramlst, edges, steps, rep, qx, maxevals)
    else:
        result = nm_noparams(q_func, init_param, problem_size, edges, steps, rep, qx, maxevals)

    total_time = time.time()  - start_time # end_time - start_time - (end_sim - start_sim)
    # print("With parameters: ", withparams)
    # print("Total execution time: '{0}' seconds and time taken by simulation was: '{1}' seconds".format(total_time, total_simulation_time))
    # print("Timecheck: %s" %(time.time() - start_time))
    line = (time.strftime("%d/%m/%y %H:%M:%S") + "\t" +  # +"Total time:\t"
                                            str(total_time) + "\t" +  #Compile time: \t" + 
                                            str(total_compile_time) +"\t" + # Simulation time:\t" + 
                                            str(total_simulation_time) + "\t" +  # problem_size:\t" + 
                                            str(problem_size) + "\t" +  #steps: \t" + 
                                            str(steps) + "\t" + # maxiter:\t" + 
                                            str(maxevals) + "\t" +  #nit:\t" + 
                                            str(result[3]) + "\t" +  #nfev:\t" + 
                                            str(result[4]) + "\t" +  #params?\t + 
                                            str(withparams) + "\n" )
    # open a file to append
    with open(os.path.join(curdir, filename), 'a') as f:
        f.write(line)
    print(line)