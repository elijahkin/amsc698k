#######################################################################
#
# QAOA helper functions for pennylane 
#   (see qaoa_utilities_intro.ipynb for examples)
#  
#######################################################################
# standard cost and mixer Hamiltonians, QAOA circuit definition, optimizer loop
#    for pennylane framework
#
# qml_standard_cost_h(weight_dict,     # dict of weights of form { (qa,qb): coeff} or { (qa): coeff } 
#                        **kwargs)     # currently not used
#      returns cost Hamiltonian based on weight_dict of form { (qa,qb): coeff} or { (qa): coeff }
#        setting PauliZ(qa) or PauliZ(qa)@PauliZ(qb)
#
# qml_general_cost_h(obs_weight_dict, # dict of weights of form { ('Oqa','Oqb'): coeff} or { ('Oqa'): coeff}
#        # or {('O_aq','O_qb'): coeff} or {('O_qa'): coeff} or { (qa,qb): coeff} or { (qa): coeff } for 'O'=Z
#                              #  where 'qa,qb' are qubits and 'O' operator ('H','X','Y','Z','S','T','V'...)
#                              #      (S=sqrt(Z), T=sqrt(S), V=sqrt(X)) 
#                        **kwargs)    # currently not used (may add controlled operators?)
#      returns cost Hamiltonian as sum of coeff*O(qa) or coeff*O(qa)@O(qb)
# 
# qml_standard_mixer_h(nwires)
#      returns mixer Hamiltonian as PauliX(qa) for each qubit
#
# qml_qaoa_initial_h(nwires, topo="x+")  #<<< not worksing as is!
#
# qml_qaoa_circuit(params, cost_h, mixer_h=None, initial_state="x+", nlayers=1)
#     note: calling 'qml_qaoa_initial_h' doesn't work yet, so 'qml_qaoa_circuit'
#           also initializes the states (same as 'topo' argument in 'qml_qaoa_initial_h')
#           initial_state="x+" or "+"  initialize all qubits to |+> state  (default)
#           initial_state="-" or "x-"  initialize all qubits to |-> state
#           initial_state="W"          initialize to superposion of W-states
#                                      (W = |100..> + |010..> + |001..> + ...)
#           initial_state="rand"       initialize qubits with random rotation about Y-axis
#           initial_state=[..,..,..]   list with rotation angle for each qubit
#
# qml_qaoa_optimize(cost_function,         # qml cost function circuit,
#                   optimizer,             # classical optimizer for this task
#                   initparams_or_nlayers, # set of initial parameters [[gammas],[alphas]] or
#                                          # number of layers (then all init. parameters are set to 0.5*np.pi
#                   max_iter=100,    # max. number of iterations (unless |diff_costs|<epsilon)
#                   epsilon=0.0,     # early termination if |diff_costs|<epsilon (epsilon=0 forces all iterations)
#                   printcycle=10,   # printout every N iterations
#                   outdata=1)       # =1: returns last parameter set and list of costs for each iteration
#                                    # =2: returns lists of all parameter sets and costs for each iteration
#
########################################

import math
import numpy as np
from datetime import datetime, timedelta
import pennylane as qml
from pennylane import qaoa
from pennylane import numpy as pnp

########

def qml_standard_cost_h(weight_dict, **kwargs):
    # skiplist: discard offsets
    if "skiplist" in kwargs.keys():
        skiplist = kwargs["skiplist"]
    else:
        skiplist = ["offset", "constant"]
    coeffs = []
    obs = []
    for k,w in weight_dict.items():
        if isinstance(k, str) and k in skiplist:
            continue
        coeffs.append(w)
        if len(k) == 1:
            obs.append(qml.PauliZ(wires=k[0]))
        else:
            obs.append(qml.PauliZ(wires=k[0])@qml.PauliZ(wires=k[1]))
    return qml.Hamiltonian(coeffs, obs)

########

def qml_general_cost_h(obs_weight_dict, **kwargs):
    
    # allowed operators (including inverse (2nd char is 'i' like 'Si')) 
    obslist = ['H','h','Z','z','Y','y','X','x','S','T','V']
    # gate dict 
    obsdict = {'Z': qml.PauliZ, 'Y': qml.PauliY, 'X': qml.PauliX, 
               'H': qml.Hadamard, 'S': qml.S, 'T': qml.T, 'V': qml.SX }
    # skiplist: discard offsets
    if "skiplist" in kwargs.keys():
        skiplist = kwargs["skiplist"]
    else:
        skiplist = ["offset", "constant"] 
    # epsilon: discard small coefficients to reduce the number of gates 
    epsilon = kwargs["epsilon"] if "epsilon" in kwargs.keys() else 0.0

    # lists of operators and coefficients as input to qml.hamiltonian()
    obs = []
    coeffs = []
    sep = None
    for k,w in obs_weight_dict.items():
        if abs(w) < epsilon:
            continue
        lk = k if isinstance(k, (list,tuple)) else [k]
        if str(lk[0]) in skiplist:
            continue
        if len(lk) == 1 and str(lk[0]).find("*") > 0:
            a = lk[0].split("*")
        else:
            a = lk
        if sep is None and isinstance(a[0], str):
            sep = "_" if a[0].find("_") > 0 else ""
        qA, qB = -1,-1 
        if isinstance(a[0], int) or ( isinstance(a[0], str) and a[0][0].isdigit() ):
            obsA = 'Z'
            qA = int(a[0]) 
        elif a[0][0] in obslist:
            b = a[0].split(sep,1)
            obsA = b[0].upper()
            if b[1][0].isdigit():
                qA = int(b[1])
            else:
                qA = int(b[1][1:])
                obsA += b[1][0]
        if len(a) == 1:
            coeffs.append( w )
            if len(obsA) == 1:
                obs.append( obsdict[obsA](wires=qA) )
            else:
                obs.append( Adjoint(obsdict[obsA[0]](wires=qA)) )
        else:
            if isinstance(a[1], int) or ( isinstance(a[1], str) and a[1][0].isdigit() ):
                obsB = 'Z'
                qB = int(a[1])           
            elif a[1][0] in obslist:
                b = a[1].split(sep,1)
                obsB = b[0].upper()
                if b[1][0].isdigit():
                    qB = int(b[1])
                else:
                    qB = int(b[1][1:])
                    obsB += b[1][0]
            coeffs.append( w )
            obs.append( obsdict[obsA](wires=qA) @ obsdict[obsB](wires=qB) )
    return qml.Hamiltonian(coeffs, obs)

########

def qml_standard_mixer_h(nwires):
    coeffs = [1]*nwires
    obs = [qml.X(wires=j) for j in range(nwires)]
    return qml.Hamiltonian(coeffs, obs)

########

def qml_qaoa_initial_h(nwires, topo="x+"):
    """
  qml_qaoa_initial_h(nwires, topo)
    returns circuit for state initialization
    input: nwires:  number of qubits to initialize            (int, required)
           topo:    topology of initialization     (string or list, optional)
              topo="+" or "x+"  initialize all qubits to |+> state  (default)
              topo="-" or "x-"  initialize all qubits to |-> state
              topo="W"     initialize to superposion of W-states
                            (W = |100..> + |010..> + |001..> + ...)
              topo="rand"  initialize qubits with random rotation about Y-axis
              topo=number  initialize all qubits with same rotation angle about Y-axis
              topo=[...]   list with rotation angles for each qubit
    """
    vals = []
    if not isinstance(topo, str):
        if isinstance(topo, dict):
            vals = list(topo.values())
        elif isinstance(topo, (int,float)):
            vals = [topo] * nwires
        else:
            vals = topo
        topo = "rand"

    if topo == "W" or topo == "w":
        for j in range(nwires-1,0,-1):
            th = -2. * math.asin(1./math.sqrt(j+1))
            qml.CRY(th, wires=[j, j-1])
            qml.CNOT(wires=[j-1, j])

    elif topo[0] in ["X","x","+","-"]:
        if len(topo) > 1:  topo[0] = topo[1]
        for j in range(nwires):
            if topo[0] == "-":
                qml.X(wires=j)
            qml.Hadamard(wires=j)

    elif topo == "rand":
        if len(vals) == 0:
            vals = np.random.uniform(low=0.5*np.pi, high=2.5*np.pi, size=nwires)
        for j in range(nwires):
            qml.RY(vals[j], wires=j)
#???
    return qml

#######

def qml_qaoa_circuit(params, cost_h, mixer_h=None, initial_state="x+", nlayers=1):

    if type(cost_h) == "function":
       my_cost_h = cost_h()
    else:
       my_cost_h = cost_h
    nwires = len(my_cost_h.wires)
    #print("C", my_cost_h)

    if mixer_h is None:
        my_mixer_h = qml_standard_mixer_h(nwires)
    elif type(mixer_h) == "function":
        my_mixer_h = mixer_h()
    else:
        my_mixer_h = mixer_h
    #print("M", my_mixer_h)

    coeffs = []
    if not isinstance(initial_state, str):
        coeffs = list(initial_state.values()) if isinstance(initial_state, dict) else initial_state
        initial_state = "rand"

    if initial_state == "W" or initial_state == "w":
        for j in range(nwires-1, 0, -1):
            th = -2. * math.asin(1./math.sqrt(j+1))
            qml.CRY(th, wires=[j, j-1])
            qml.CNOT(wires=[j-1, j])
    elif initial_state[0] in ["X","x","+","-"]:
        if len(initial_state) > 1:  initial_state = initial_state[1:]
        for j in range(nwires):
            if initial_state[0] == "-":
                qml.X(wires=j)
            qml.Hadamard(wires=j)
    elif initial_state[0] == "r":
        if len(coeffs) == 0:
            coeffs = np.random.uniform(low=0.5*np.pi, high=1.5*np.pi, size=nwires)
        for j in range(nwires):
            qml.RY(coeffs[j], wires=j)

    def qaoa_layer(gamma, beta):
        qaoa.cost_layer(gamma, my_cost_h)
        qaoa.mixer_layer(beta, my_mixer_h)

    qml.layer(qaoa_layer, nlayers, params[0], params[1])

########

def qml_qaoa_optimize(cost_function, optimizer, initparams_or_nlayers,
                      max_iter=100, epsilon=0.0, printcycle=10, outdata=1, **kwargs):

    #stop optimization if cost doesn't change more than epsilon over 3 iterations
    
    # initialize parameters
    if isinstance(initparams_or_nlayers, int):
        nlayers = initparams_or_nlayers
        initparams = np.repeat(0.5*np.pi, 2*nlayers).reshape(2,nlayers).tolist()
    elif isinstance(initparams_or_nlayers, np.ndarray):
        initparams = initparams_or_nlayers.tolist()
    elif isinstance(initparams_or_nlayers, dict):
        initparams = list(initparams_or_layers.values())
    else:
        initparams = initparams_or_nlayers
    if not isinstance(initparams[0], list):
        nlayers = len(initparams)//2
        initparams = [ initparams[:nlayers] , initparams[nlayers:] ]
    params = pnp.array(initparams, requires_grad=True)
    if printcycle > 0:
        print("Init params:", np.round(initparams,4).tolist())
        print(f"  Max. Iterations={max_iter}  unless |diff(costs)|<{epsilon}")
        drawcircuit = outdata<0

    outdata = abs(outdata)
    mycosts = [6000, 5000, 4000]
    costs = [5000,]
    if outdata > 1:
        allparams = [ initparams ]

    if printcycle > 0:
        start_time = datetime.now()
        print("\nStart of optimization loop:", start_time)

    for i in range(1, max_iter+1):

        params, prev_cost = optimizer.step_and_cost(cost_function, params, **kwargs)

        if outdata > 1:
            allparams.append( np.round(params,4).tolist() )
        costs.append(np.round(prev_cost,4))
        mycosts[i%3] = abs(costs[i-1] - prev_cost)

        if printcycle > 0 and i%printcycle == 0:
            print(f"Iteration {i}: cost={costs[-1]} (cost diff={round(mycosts[i%3],4)})")
            if drawcircuit:
                print(qml.draw(cost_function, expansion_strategy="device")(params))

        if mycosts[0]<epsilon and mycosts[1]<epsilon and mycosts[2]<epsilon:
            if printcycle > 0:
                print("Early termination at iteration",i)
            break

    # add costs after last iterations
    costs.append( pnp.round(cost_function(params),4).tolist() )
    costs = costs[1:]
    # display optimized parameters and final costs
    params = np.round(params,4).tolist()
    if printcycle > 0:
        end_time = datetime.now()
        print("End of optimization loop: {} (duration: {})\n".format(end_time, end_time-start_time))
        print("Optimized parameters:", params)
        print("Minimum energy:", costs[-1])
    if outdata == 0:
        return params
    elif outdata == 1:
        return params, costs
    elif outdata == 2:
        return allparams, costs

#########################################3

