##################################################################
#
# helper functions for QAOA and Graph manipulation and plotting
#
##################################################################
# helper functions to extract information from graph objects or position or weight arrays or lists
##########
#
## for TSP: transforms between point (node) number and binary representation
##    e.g. [0,1,2] <-> '100010001' or [0,1,2] <-> [[1,0,0],[0,1,0],[0,0,1]]
#
# bitstring_to_points(bitstring,             # (0's or 1's), can be int (decimal of bitstring) or string
#                     starting_point=None)   # fixed starting_point (NOT part of the bitstring!)
#    transformation from bitstring representation to point numbers (0...N)
#
# points_to_bitstring(points,                # list (or tuple or array) of point numbers (0...N)
#                     starting_point=None,   # fixed starting_point (NOT listed in 'points'!)
#                     tsp_matrix=False)      # flag to define the output format:
#                  tsp_matrix=0 (=False): single bitstring containing the order of all points
#                  tsp_matrix=1 (=True):  list of binary lists with position of point 0 (1st row), then
#                                         position of point 1 (2nd row), etc
#                  tsp_matrix=2:  same as tsp_matrix=1 but returns numpy array
#
# get_position_dict(P_or_G=None,    # graph object or array or list or dict with node positions
#                   nodes=None,     # number of nodes or dict or list with nodenames
#                   seed=0,         # seed to create random positions if only 'nodes' is given
#                                   # skip creating random positions if seed<0; seed=0: don't set seed
#                   scale=1.0,         # scale factor for random numbers for positions
#                   rounding_digits=4, # rounding of positions
#                   verbose=1):        # optional printout
#       if possible get nodenames & positions from graph or array; otherwise create positions
#       returns node_position dict of form { node: [xposition, yposition] }
#
# get_ids_from_keys(keylist,        # list of keys (str or int)
#                   skiplist=[]):   # list of keys to skip from results (e.g. "offset")
#       extract node ids from keys (optionally discard keys (skiplist))
#       returns keydict of form { key:(nodeA,nodeB) } and namedict of form { nodeA:number }
#
# get_weight_matrix(W_or_G=None,    # graph object or array or list or dict with weights
#                   Pos=None,       # array or list or dict of positions (to calculate weights=distances)
#                   verbose=1):
#       returns weight matrix (Ising matrix) as numpy.ndarray
#
# get_normalized_weights(W_or_G,           # graph object or dict,list,array with weights
#                       scale=2*np.pi):    # scalefactor if scale<0 or normalize to scale value
#       returns dict of form {(nodeA,nodeB): coefficient} or {(nodeA): coefficient)}
#       to use directly for cost hamiltonian
#
# get_most_probable_states(probabilities,     # dictionary or list with probabilities
#                          nresults=1,        # number of most probable states
#                          minprob=0.0,       # return results for states with probability>minprob (discards 'nresults')
#                          rounding_digits=4) # rounding of probabilities for output
#       returns dict of form { state: probability } of 'nresults' most probable states
#
#
##################################################################
# helper functions for graph creation and plotting
##########
#
# create_graph(nodes=None,     # number of nodes (integer) or list of nodenames or dict
#              positions=None, # 2-dim array with x,y position per node
#              weights=None,   # weight matrix or dict
#              edges=None,     # same as weight
#              seed=0,         # seed for random position layout (if seed!=0)
#              scale=10,       # scale factor for random position layout via np.random.rand() (for seed>0)
#              verbose=False)  # verbose=1: print graph info, verbose=2: print info about array creation as well
#
# plot_graph(J_or_G,        # input: Ising matrix or Graph object
#            colors=None,   # optional: node coloring
#            cmap=None,     # optional: provide color map (e.g. "viridis", "magma", "plasma", ...)
#            plot_weights=False,  #optional: plot weights of edges (values and line width)
#            pos=None,      # optional: provide positions of nodes (otherwise: nx.spring_layout())
#            seed=0,        # optional: set random seed for nx.spring_layout()
#            filename=None, # filename to save the figure
#            **kwargs)      # additional plotting parameters
#
# plot_histogram(data,                # single dataset or list or array of datasets
#                label=None,          # label for each dataset
#                xvals=None,          # list or array of x-values (same length as datasets)
#                xlabel=None,         # label for x-axis
#                ylabel=None,         # label for y-axis
#                legend_loc="best",   # location of legend box: "best","upper right","upper left", etc)
#                filename=None,       # filename to save the figure
#                **kwargs):           # additional plot parameters (e.g. title=.., figsize=(width,height), etc)
#       simple plot function e.g. plot_histogram(costs,"Costs vs Iteration")
#       can plot multiple graphs and graphs in seperate subplots
#
# plot_probabilities(probabilities,    # dict with key=state or list of probabilities for all possible states
#                    all_states=False, # all possible states on x-axis if True
#                    filename=None,    # save figure in file
#                    **kwargs)         # additional plotting parameters
#
#
##################################################################
#### helper functions for Traveling Salesman Problem
#
# generate coefficients for TSP starting at location 0 using (N-1)**2 qubits
#      (in matrix form: rows=locations, columns=time steps; ideally initialized with superposition of W-states)
#  (based on "Important Quantum Gates for Quantum Algorithms of Travelling Salesman Problem",
#    by T.J.H.Sinaga et al, in Proceedings of ICoABCD 2023)
#
# tsp_generate_coeffs(J_or_G,           # Ising matrix (distance matrix) or Graph
#                     scale=2*np.pi,    # scale factor if scale<0 or normalize to scale value
#                     epsilon=0.0,      # set coefficient to zero if scaled coefficient < epsilon
#                     **kwargs)         # other parameters (not used yet)
#
#
## transforms between point (node) number and binary representation
##    e.g. [0,1,2] <-> '100010001' or [0,1,2] <-> [[1,0,0],[0,1,0],[0,0,1]]
#
# bitstring_to_points(bitstring,             # (0's or 1's), can be int (decimal of bitstring) or string
#                     starting_point=None)   # fixed starting_point (NOT part of the bitstring!)
#    transformation from bitstring representation to point numbers (0...N)
#      output: list of points
#
# points_to_bitstring(points,                # list (or tuple or array) of point numbers (0...N)
#                     starting_point=None,   # fixed starting_point (NOT listed in 'points'!)
#                     tsp_matrix=False)      # flag to define the output format:
#    transformation from point numbers (0...N) to bitstring (tsp_matrix=False) or matrix of binaries (tsp_matrix>0)    
#      output for tsp_matrix=0 (=False): single bitstring containing the order of all points
#                 tsp_matrix=1 (=True):  list of binary lists with position of point 0 (1st row), then
#                                         position of point 1 (2nd row), etc
#                 tsp_matrix=2:  same as tsp_matrix=1 but returns numpy array
# 
# tsp_get_most_probable_paths(probabilites,      # dict with key=state or list of probabilities for all possible states
#                             nresults=3,        # number of N best results to report in output
#                             tsp_fidelity=2,    # output tsp_matrix 
#                             list_states=False, # list of all states contributing to this output (see below)
#                             rounding_digits=4) # rounding of probabilities for output
#     returns sorted dict of results:
#       for tsp_fidelity=0: all results independent of whether they represent a 'good' path
#       for tsp_fidelity=1: results where each place was visited (even if 2 places are visited at the same time)
#       for tsp_fidelity=2: 'good' results (which contain all places (nodes) visited at different times).
#     output:
#       if list_states==False (default):
#          dict of form {'good_path': [colorlist, probability]}
#       if list_states==True:
#          dict of form {'good_path': [colorlist, probability, list_of_states]}
#        (note: good path is e.g. ['1000','0001','0100','0010'] but also ['1000','0101','0110','0001'] (where
#              first all unique timestamps are counted, the remaining are split between remaining
#              places (nodes) if they have a '1' for these timestamps: (here: ['1000','0100','0010','0001']))
#
##################################################################
##################################################################

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, OrderedDict
import math

# helper functions: 
#
## TSP transformation: points <-> bitstring representation 

def bitstring_to_points(bitstring, starting_point=None):
    """
 bitstring_to_points(bitstring, starting_point)
    transformation from bitstring representation to point numbers (0...N)
        bitstring representation (0's or 1's), this can be int (decimal of bitstring) or string  
                        or list of ints, strings or list of list of 0,1
        starting_point: number of the point if the starting_point is fixed and is NOT part of 
                        the bitstring  (default: None)
    """
    mybitstrings = []
    if isinstance(bitstring, np.ndarray):
        bitstring = bitstring.tolist()
    if isinstance(bitstring, (list,tuple)):  #list of bits or nested list (matrix)
        if isinstance(bitstring[0], (list,tuple)):
            for i in range(len(bitstring)):
                mybitstrings.append( "".join([str(b) for b in bitstring[i]]))
        else:
            for b in bitstrings:
                if isinstance(b, int):
                    mybitstrings.append( bin(b)[2:] )
                else:
                    mybitstrings.append( b[2:] if b.startswith("0b") else b )
        npoints = len(mybitstrings)
    else:
        if isinstance(bitstring, int):
            bitstring = bin(bitstring)[2:]
        if bitstring.startswith("0b"):         
            bitstring = bitstring[2:]
        npoints = int(np.sqrt(len(bitstring)))
        mybitstrings = [ bitstring[i:i+npoints] for i in range(0,len(bitstring),npoints) ]
    points = [0] * npoints
    for i,b in enumerate(mybitstrings):
        points[b.find("1")] = i
    if starting_point is not None:
        points = [starting_point] + points
    return points

def points_to_bitstring(points, starting_point=None, tsp_matrix=False):
    """
 points_to_bitstring(points, starting_point, tsp_matrix)
    transformation from point numbers (0...N) to bitstring representation
       points:  list (or tuple or array) of point numbers (0...N)
       starting_point: number of the point if the starting_point is fixed and is NOT part of
                        the bitstring  (default: None)
       tsp_matrix:      flag to define the output format
          tsp_matrix=0 (=False): single bitstring containing the order of all points
          tsp_matrix=1 (=True):  list of binary lists with position of point 0 (1st row), then
                                 position of point 1 (2nd row), etc
          tsp_matrix=2:  same as tsp_matrix=1 but returns numpy array
    """
    if isinstance(points, np.ndarray):
        points = points.tolist()
    if not isinstance(points, (list,tuple)): 
        return None
    tsp_matrix = int(tsp_matrix)
    if starting_point is not None:
        points = [starting_point] + points
    npoints = len(points)
    mat = np.zeros((npoints,npoints), dtype=int)
    for i,p in enumerate(points):
        mat[p][i] = 1
    if tsp_matrix > 0:
        return mat.tolist() if tsp_matrix==1 else mat
    return "".join([ str(mat[i][j]) for i in range(mat.shape[0])
                    for j in range(mat.shape[1]) ])

###########

def get_position_dict(P_or_G=None,    # graph object or array or list or dict with node positions
                      nodes=None,     # number of nodes or dict or list with nodenames
                      seed=0,         # seed to create random positions if only 'nodes' is given
                                      # skip creating random positions if seed<0; no seed set if seed=0
                      scale=1.0,          # scale factor for random numbers for positions
                      rounding_digits=3,  # rounding of x,y position values 
                      verbose=1):         # optional printout

    # (a) P_or_G is graph object or dict or array or list of positions (or list of [nodename,position])
    pos = {}
    if P_or_G is None and nodes is not None:
        P_or_G = nodes
        nodes = None
    if P_or_G is not None:
        if isinstance(P_or_G, nx.classes.graph.Graph):
            pos = nx.get_node_attributes(P_or_G, "pos")
        elif isinstance(P_or_G, dict):
            pos0 = list(P_or_G.values())[0]
            if isinstance(pos0, (list,tuple)):    # otherwise dict with nodenames?
                pos = P_or_G
            else:
                nodes = P_or_G
        elif isinstance(P_or_G, np.ndarray):
            plist = P_or_G.tolist() 
            if nodes is not None and isinstance(nodes, (list,tuple)):
                pos = { n:p for n,p in zip(nodes, plist) }
            else:
                pos = { j:p for j,p in enumerate(plist) }
        elif isinstance(P_or_G, list):
            if isinstance(P_or_G[0], (list,tuple)):
                if len(P_or_G[0]) == 3:
                    for p in P_or_G:
                        pos[p[0]] = ( p[1], p[2] )
                else:
                    if isinstance(P_or_G[0][1], (list,tuple)):
                        for p in P_or_G:
                            pos[p[0]] = p[1]
                    else:
                        pos = { j:p for j,p in enumerate(P_or_G) }
            else:
                # otherwise list of nodenames?
                nodes = P_or_G
        elif isinstance(P_or_G, int):
            nodes = P_or_G                   # assume number of nodes provided

    if len(pos) > 0:
        return { j:np.round(p,rounding_digits).tolist() for j,p in pos.items() }

    if nodes is None or seed < 0:
        if verbose > 0:
            if seed < 0 and nodes is not None:
                print("*** get_position_dict: Creation of random position skipped as seed<0")
            else:
                print("*** get_position_dict: Provide Graph object or list or dict or array with positions")
        return None

    # (b) if position not provided, create random positions for numNodes
    if verbose > 1:
        print(f"* get_position_dict: Create random positions for {nodes} nodes")
    if isinstance(nodes, dict):
        nodelist = list(nodes.values())
        if isinstance(nodelist[0], int) and sorted(nodelist)[-1] == len(nodes)-1:
            nodelist = list(nodes.keys())
    elif isinstance(nodes, int):
        nodelist = list(range(nodes))
    elif isinstance(nodes, np.ndarray):
        nodelist = nodes.tolist()
    else:
        nodelist = nodes 
    if seed > 0:
        np.random.seed(seed)
    pos = np.round( np.random.rand(len(nodelist),2) * scale, rounding_digits)
    return { j:p for j,p in zip(nodelist,pos.tolist()) }

##########################

def get_ids_from_keys(keylist,       # list of keys (str or int)
                      skiplist=[]):  # list of keys to skip from results (e.g. "offset")
    """
 get_ids_from_keys(keylist, skiplist) 
   extract node ids from keys (optionally discard keys (skiplist))
   returns keydict of form { key:(nodeA,nodeB) } and namedict of form { nodeA:number }  
    """
    if isinstance(keylist, dict):
        keylist = list(keylist.keys())
    idx = []
    keydict = {}
    sep0 = None; sep1 = None
    for k in keylist:
        if isinstance(k, (list,tuple)):
            if len(k) == 1:
                idx.append( k[0] )
                keydict[k] = ( k[0], )
            else:
               idx.extend([ k[0], k[1] ])
               keydict[k] = ( k[0], k[1] )
        elif isinstance(k, int) or ( isinstance(k, str) and k.isdigit() ):
            idx.append( k )
            keydict[k] = ( int(k), )
        else:
            if k in skiplist:
                keydict[k] = ()
                continue
            if sep0 is None:
                for w in "*@&-,:":
                    if k.find(w) > 0:
                        sep0 = w
                        break
            l = k.split(sep0)
            m = l if isinstance(l, str) else l[0]
            llen = 1 if isinstance(l, str) else len(l)
            if sep1 is None:
                if m.isdigit():
                     sep1 = "?"
                else:
                    for w in "_xyzh^-":
                        if m.find(w) > 0:
                            sep1 = w
                            break
            if llen == 1:
                idx.append( m.split(sep1)[-1] )
                keydict[k] =  ( m.split(sep1)[-1], ) 
            else:
                m = l[0].split(sep1)[-1]; n = l[1].split(sep1)[-1]
                idx.extend([ m, n ])
                keydict[k] = ( m, n )

    myids = list(dict.fromkeys(idx).keys())
    # special cases: myids are digits or strings containing a sequence of ints
    #  e.g.  myids=['1','4','3','2'] -> return sorted list
    if isinstance(myids[0], int):
        namedict = { k:j for j,k in enumerate(sorted(myids)) }
    else:
        tmp = [ int(m) for m in myids if isinstance(m, str) and m.isdigit() ]
        stmp = sorted(tmp)
        if len(tmp) == len(myids) and stmp[-1] == stmp[0] + len(tmp) -1:
            namedict = { str(k):j for j,k in enumerate(stmp) }
        else:
            namedict = { k:j for j,k in enumerate(myids) }

    return keydict, namedict

##############################

def get_weight_matrix(W_or_G=None, Pos=None, verbose=1):

    # (a) position matrix or dict or list of coordinates given
    if W_or_G is None:
        if Pos is None:
            if verbose > 0:
                print("Provide Graph object or list or dict or array with edge weights")
            return None
        node_dict = get_position_dict(Pos)
        numNodes = len(node_dict)
        ids = { k:j for j,k in enumerate(list(node_dict.keys())) }
        W = np.zeros( (numNodes,numNodes) )
        for nA,pA in node_dict.items():
            for nB,pB in node_dict.items():
                if nA != nB:
                    W[ids[nA]][ids[nB]] = np.sqrt((pA[0] - pB[0])**2 + (pA[1] - pB[1])**2)
        return W

    # (b) full weight matrix
    if isinstance(W_or_G, nx.classes.graph.Graph):
        return nx.to_numpy_array(W_or_G)
    if isinstance(W_or_G, np.ndarray):
        if W_or_G.shape[0] == W_or_G.shape[1]:
            return W_or_G
        else:
            if verbose > 0:
                print("Unexpected size of weight matrix! Use 'Pos' input argument?")
            return None
    if isinstance(W_or_G, list) and isinstance(W_or_G[0], list) and len(W_or_G[0]) == len(W_or_G):
        return np.array(W_or_G)

    # (c) list or dict with weights
    mykeys = []
    myweights = []
    if isinstance(W_or_G, dict):
        mykeys = list(W_or_G.keys())
        myweights = list(W_or_G.values())
    elif isinstance(W_or_G, list):
        if isinstance(W_or_G[0], (list, tuple)): 
            if ( len(W_or_G[0]) == 2 and  
                (( isinstance(W_or_G[0][0], str) and isinstance(W_or_G[0][1], (int,float)) ) or
                 ( isinstance(W_or_G[0][0], int) and isinstance(W_or_G[0][1], float) )) ):
                for entry in W_or_G:
                     mykeys.append( entry[0] )
                     myweights.append( entry[1] )
            else:
                for entry in W_or_G:
                     mykeys.append( (entry[0],entry[1]) )
                     myweights.append( 1.0 if len(entry)==2 else entry[2] )
        elif instance(W_or_G[0], str):
            for entry in W_or_G:
                 mykeys.append( entry )
                 myweights.append( 1.0 )
    if len(mykeys) == 0:
        if verbose > 0:
            print("Cannot find node IDs in list")
        return None

    # extract nodes from keys
    skiplist = []
    for k in mykeys:
        m = str(k) if isinstance(k, (int,str)) else str(k[0])
        if m.startswith("off") or m.startswith("const"):
           skiplist.append(k)
    mydict, ids = get_ids_from_keys(mykeys, skiplist)

    numNodes = len(ids)
    W = np.zeros((numNodes,numNodes))
    for k,w in zip(mykeys, myweights):
        n = mydict[k]
        if len(n) > 1:
            W[ids[n[0]]][ids[n[1]]] = w
            W[ids[n[1]]][ids[n[0]]] = w
        elif len(n) == 1:
            W[ids[n[0]]][ids[n[0]]] = w    
    return W

###########################################

def get_normalized_weights(W_or_G,           # graph object or dict,list,array with weights
                          scale=2*np.pi):    # scalefactor if scale<0 or normalize to scale value 

# output: dict of form { (nodeA,nodeB): coefficient } or { (nodeA): coefficient) }

    weight_dict = {}
    # (a) graph object: collect weights of nodes and edges
    if isinstance(W_or_G, nx.classes.graph.Graph):
        weights = list(nx.get_edge_attributes(W_or_G, "weight").values())
        ids = { k:j for j,k in enumerate(list(W_or_G.nodes())) }
        posweights = list(nx.get_node_attributes(W_or_G, "weight").values())
        weights.extend(posweights)
        if len(weights) == 0:
            for i,j in W_or_G.edges():
                weight_dict[ ( ids[i], ids[j] ) ] = 1.0
        else:
            scale = np.max(weights)/scale if scale > 0 else abs(scale)
            if len(posweights) > 0:
                for i,w in nx.get_node_attributes(W_or_G, "weight").items():
                    weight_dict[ ( ids[i], ) ] = w/scale
            for k,w in nx.get_edge_attributes(W_or_G, "weight").items():
                a,b = k
                weight_dict[ ( ids[a], ids[b] ) ] = w/scale
        return weight_dict

    # (b) weight_matrix (array or nested lists) -> mylen=number_of_rows
    if isinstance(W_or_G, np.ndarray) and W_or_G.shape[0] == W_or_G.shape[1]:
        mylen = W_or_G.shape[0]
    elif isinstance(W_or_G, (list,tuple)) and isinstance(W_or_G[0], (list,tuple)):
        mylen = len(W_or_G) if len(W_or_G)==len(W_or_G[0]) else 1
    elif isinstance(W_or_G, dict):
        mylen = 0
    if mylen > 1:
        scale = np.max(W_or_G)/scale if scale > 0 else abs(scale)
        for i in range(mylen):
            if W_or_G[i][i] != 0.0:
                weight_dict[ ( i, ) ] = W_or_G[i][i]/scale
            for j in range(i+1, mylen):
                if W_or_G[i][j] != 0.0:
                    weight_dict[ ( i, j ) ] = W_or_G[i][j]/scale
        return weight_dict

    # (c) dict or list of form [[nodeA_nodeB,coeff],..] or [[nodeA,nodeB,coeff],..] or [[nodeA,coeff],..]
    # remove entries "offset" and "constant"
    mykeys = []; weights = []
    if mylen==0:
        for k,w in W_or_G.items():
            m = str(k) if isinstance(k, (str,int)) else str(k[0])
            if m.startswith("off") or m.startswith("const"):
                continue
            weights.append( w )
            mykeys.append( k )
    else:
        for entry in W_or_G:
            m = str(entry[0])
            if m.startswith("off") or m.startswith("const"):
                continue
            weights.append( entry[-1] )
            mykeys.append( entry[:-1] )
    scale = np.max(weights)/scale if scale > 0 else abs(scale)
    
    iddict, ids = get_ids_from_keys(mykeys)    
    for k,w in zip(mykeys,weights):
        m = iddict[k]
        if len(m) == 1:
            weight_dict[ ( ids[m[0]], ) ] = w/scale
        elif len(m) == 2:
            weight_dict[ ( ids[m[0]], ids[m[1]] ) ] = w/scale
    return weight_dict 

###########################################

def get_most_probable_states(probabilities,   # dictionary or list with probabilities
                             nresults=1,      # number of most probable states
                             minprob=0.0,     # return results for states with probability>minprob (discards 'nresults')
                             rounding_digits=4, # rounding of results
                             **kwargs):
    """
 get_most_probable_states( probabilities, nresults, minprob, rounding_digits)
   returns dict {state: probability} for 'nresults' most probable states 
   input: probabilities:  dictionary with state as key and prob. as value;
                          or list with probs for all possible states
          nresults:       number of most probable states for output                    (default: 5)
          minprob:        return results for states with probability>minprob (discards 'nresults'; default: 0.0)
          rounding_digits: rounding of results       (default: 4)
     """

    if len(kwargs) > 0 and "tsp_matrix" in kwargs.keys():
        print("For TSP related results please use 'tsp_get_most_probable_paths()'")

    if len(probabilities) < nresults:
        print("Requested number of most probable results ({}) > elements in probability array ({})".format(
               nresults,len(probabilities)))
        nresults = len(probabilities)
    if minprob > 0:
        minprob += 0.000001
        nresults = len(probabilities)

    mylocmax = {}
    mydict = None
    if isinstance(probabilities, dict):
        mydict = probabilities
    else:
        # list or array of samples?
        samples = []
        tmp = probabilities.tolist() if isinstance(probabilities, np.ndarray) else probabilities
        if isinstance(tmp[0], (list,tuple)):
             for vals in tmp:
                 samples.append("".join([ str(i) for i in vals ]))
        elif isinstance(tmp[0], str):
             for vals in tmp:
                 samples.append( vals )
        if len(samples) > 0:
            mydict = defaultdict(int)
            for entry in samples:
                 mydict[entry] += 1

    if mydict is not None:
        keylen = 1
        for j in range(nresults):
            locmax, probmax = 0, -1.0
            for k,v in mydict.items():
                if j == 0:
                    ck = bin(k)[2:] if isinstance(k, int) else k
                    if len(ck) > keylen:
                        keylen = len(ck)
                if v > probmax:
                    locmax, probmax = k, v
            if probmax > minprob:
                mylocmax[locmax] = probmax
                mydict[locmax] = -1.0
        for k,v in mylocmax.items():
            mydict[k] = v

    else:
        if len(tmp)>0:
            keylen = math.ceil(math.log2(len(tmp)))
        for _ in range(nresults):
            probmax = float(np.max(probabilities))
            locmax  = int(np.argmax(probabilities))
            if probmax > minprob:
                mylocmax[locmax] = probmax
                probabilities[locmax] = -1.0
        for k,v in mylocmax.items():
            probabilities[k] = v

    locmaxarr = {}
    for k,v in mylocmax.items():
        if isinstance(k, int):
            locmaxarr[ bin(k)[2:].zfill(keylen) ] = round(v, rounding_digits)
        else:
            locmaxarr[ k.zfill(keylen) ] = round(v, rounding_digits)

    return locmaxarr


##################################################################

def create_graph(nodes=None,     # number of nodes (integer) or list of nodenames or dict
                 positions=None, # 2-dim array with x,y position per node
                 weights=None,   # weight matrix or dict
                 edges=None,     # same as weight
                 seed=0,         # seed for random position layout (if seed>0); nx.spring_layout() if seed<0
                 scale=10.,      # scale factor for random position layout via np.random.rand() (for seed>0)
                 verbose=False,
                ):
    """
  create_graph( nodes, positions, weights, edges, seed, scale, verbose )
    create networkx graph object for various types of input data 
              (at least one of 'nodes', 'positions', 'weights' or 'edges' must be provided)
    input: nodes:      number of nodes or list or dict of nodenames    (int or list or dict, optional)
           positions:  2-dim array with x,y position per node (array or nested list or dict, optional)
           weights:    weight matrix or dict with (weighted) edges  (2-dim array or nested list or dict, optional)
           edges:      same as weight
           seed:       seed for random node position (if seed>0) or nx.spring_layout (if seed<0)   (int, optional)
           scale:      multiplier for random position layout via np.random.rand()           (default=10, optional)
           verbose:    info printout (=True (or 1): graph creation; =2: also generation of missing arrays)
                                                                                           (bool or int, optional)
    output: Graph object
    """
    verbose = int(verbose)

    # nodenames and positions
    node_dict = get_position_dict(P_or_G=positions, nodes=nodes, seed=seed, scale=scale, verbose=verbose)

    # weight matrix (edges)
    if edges is not None:
        if weights is None:
            weights = edges
        else:
            if isinstance(weights, np.ndarray):
                weights = weights.tolist()
            if len(weights) == len(edges):
                myweights = {}
                for key,w in zip(edges, weights):
                    myweights[key] = w
                weights = myweights
        if verbose > 1:
            print("edges - weights?:", weights)

    weight_matrix = get_weight_matrix(W_or_G=weights, Pos=node_dict, verbose=verbose)
    
    if node_dict is not None:
        numNodes = len(node_dict)
        nodes = list(node_dict.keys())
    elif weight_matrix is not None:
        numNodes = weight_matrix.shape[0]
        mykeys = []
        if isinstance(weights, dict):
           mykeys = list(weights.keys())
        elif isinstance(weights, list) and isinstance(weights[0], (list,tuple)):
            if len(weights[0]) <=3:           # last value is weight
                for k in weights:
                    mykeys.append( k[:-1] )
        if len(mykeys) > 0:
           keydict,ids = get_ids_from_keys(mykeys)
           nodes = list(ids.keys())
    elif nodes is not None:
        numNodes = nodes if isinstance(nodes, int) else len(nodes)

    if numNodes == 0:
        print("*** NO GRAPH CREATED: neither nodes nor positions nor weights given ***")
        return None

    if nodes is None or len(nodes) == 1:
        nodes = list(range(numNodes))

    gr = nx.Graph()

    if verbose:
        print("graph -nodes:",nodes)
        print("graph -node_dict:",node_dict)
        print("graph -weight_matrix:",weight_matrix)

    if node_dict is None: 
        for i in nodes:
            gr.add_node(i)
        node_dict = nx.spring_layout(gr, seed=abs(seed))
        if verbose:
            print("graph -seed<=0: node_dict",node_dict)
    else:
        for i,pos in node_dict.items():
            gr.add_node(i, pos=pos)
        if verbose > 1:
            print("graph -node_dict!=None: pos",nx.get_node_attributes(gr,"pos"))

    if weight_matrix is None:
        skipedges = []
        for k0,p0 in node_dict.items():
            skipedges.append(k0)
            for k1,p1 in node_dict.items():
                if k1 not in skipedges:
                    gr.add_edge(k0, k1, 
                                weight=np.sqrt( (p0[0]-p1[0])**2 + (p0[1]-p1[1])**2 ))
        if verbose > 1:
            print("graph -weight=None: add_edge", nx.get_edge_attributes(gr,"weight"))
    else:
        for i,ni in enumerate(gr.nodes()):
            for j,nj in enumerate(gr.nodes()):
                if weight_matrix[i][j] != 0.0:
                    gr.add_edge(ni, nj, weight=weight_matrix[i][j])
        if verbose >1:
            print("graph -weight_matrix: add_edge", nx.get_edge_attributes(gr,"weight"))

    return gr

#######################################################################################

def plot_graph(J_or_G,        # input: Ising matrix or Graph object
               colors=None,   # optional: node coloring
               cmap=None,     # optional: provide color map (e.g. "viridis", "magma", "plasma", ...)
               plot_weights=False,  #optional: plot weights of edges (values and line width)
               pos=None,      # optional: provide positions of nodes (otherwise: nx.spring_layout())
               seed=0,        # optional: set random seed for nx.spring_layout()
               filename=None, # filename to save the figure
               **kwargs):     # additional plotting parameters
    """ 
  plot_graph( J_or_G, colors, cmap, plot_weights, pos, seed, filename, **kwargs)
    plot graph with partition (colors)
    input: J_or_G:  Ising (dist matrix) or Graph object
           colors:  node coloring                     (array or list or dict or bitstring, optional)
                     (first node = highest bit in bitstring or first entry in list or array)
           cmap:    provide color map (e.g. "viridis", "magma", "plasma", ...)    (string, optional)
           plot_weights: plot weights of edges (values and line width)           (boolean, optional)
           pos:          provide positions of nodes         (array or nested list or dict, optional)
           seed:         set random seed for nx.spring_layout() in case that pos=None (int, optional)
           filename:     filename to save the figure                              (string, optional)
           kwargs:       additional plotting parameters                     (dict entries, optional)
    output: none
    """
    if isinstance(J_or_G, nx.classes.graph.Graph):
        graph = J_or_G
        N = graph.number_of_nodes()
    else:
        # define graph
        graph = nx.Graph()
        if isinstance(J_or_G, (list,tuple)):
            J_or_G = np.array(J_or_G)
        if J_or_G.ndim != 2:
            print("*** plot_graph: weight matrix is not 2-dim array")
            return 
        N = J_or_G.shape[0]
        for ii in range(0, N):
            for jj in range(ii + 1, N):
                if J_or_G[ii][jj] != 0:
                    graph.add_edge(ii, jj, weight=J_or_G[ii][jj])

    nodeID = { k:j for j,k in enumerate(list(graph.nodes())) }
    if pos is None:
        pos = nx.get_node_attributes(graph, "pos")
        if len(pos) == 0 or seed > 0: 
            if seed > 0:
                pos = nx.spring_layout(graph, seed=seed)
            else:
                pos = nx.spring_layout(graph)
    mypos = np.zeros((N,2)).tolist()
    if isinstance(pos, dict):
        for k,p in pos.items():
            mypos[nodeID[k]] = p
    else:
        if isinstance(pos, np.ndarray):
            mypos = pos.tolist()
        else:
            mypos = pos

    # input: bitstring; output: list of colors (numbers)
    def bitstring_to_colors(bitstring, nNodes):
        if not isinstance(bitstring, str) or len(bitstring) == 0:
            return []
        if len(bitstring) == 1:
            return [bitstring] * nNodes
        if bitstring[:2] == '0b':
            bitstring = bitstring[2:]
            if len(bitstring) == 1:
                return [bitstring] * nNodes
        if bitstring[:2] not in ['00','01','10','11']:
            return [bitstring] * nNodes
        if len(bitstring) < N:
            return [ int(b) for b in list(bitstring) ] + [0] * (N-len(bitstring))
        if len(bitstring) == N*N or len(bitstring) == (N-1)*(N-1):
            mylen = N if len(bitstring)==N*N else (N-1)
            mycols = [] if mylen==N else [0,]
            for j in range(0, len(bitstring), mylen):
                 onerow = bitstring[j:j+mylen]
                 mycols.append( onerow.find('1') + N-mylen )
            return mycols
        return [ int(b) for b in list(bitstring) ]

    # define color scheme
    if colors is None:
        mycolors = nx.get_node_attributes(graph, "color")
    if isinstance(colors, dict):
        key0 = list(colors.keys())[0]    # most probable state
        mycolors = bitstring_to_colors(key0, N)
        if len(mycolors) == 0:
           mycolors = list(colors.values())
    elif isinstance(colors, (list,tuple)):
        mycolors = colors.copy()
    elif isinstance(colors, int): 
        mycolors = [ colors for _ in range(N) ]
    else:
        mycolors = bitstring_to_colors(colors, N)
    if len(mycolors) == 0:
        mycolors = [0 for _ in range(N)]

    if isinstance(mycolors[0], int):
        maxcolors = np.max(mycolors)
    else:
        maxcolors = 1000   # do not use colorlist
    if maxcolors > 7 or cmap is not None:
        nodecolors = [ mycolors[j] for j in range(N) ]
        if maxcolors < 1000 and cmap is None:
            cmap = "viridis"   #cmap="magma" or "plasma"?
    else:
        colorlist = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                     "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
        nodecolors = [colorlist[mycolors[j]] for j in range(N) ]

    if plot_weights:
        # get unique weights
        all_weights = list(nx.get_edge_attributes(graph, "weight").values())
        if len(all_weights) == 0:
            all_weights = [1.] * N
        unique_weights = list(set(all_weights))
        # plot the edges - one by one
        for weight in unique_weights:
            # form a filtered list with just the weight you want to draw
            weighted_edges = [(node1, node2)
                              for (node1, node2, attr) in graph.edges(data=True)
                              if attr["weight"] == weight ]
            # width = relative weight
            width = weight * N * 10.0 / sum(all_weights)
            nx.draw_networkx_edges(graph, pos, edgelist=weighted_edges, width=width)
        minweight = np.min(np.abs(all_weights))
        digits = 1
        while minweight * 10**digits < 10.0:
            digits +=1
        labels = {e:round(float(w),digits) for e,w in nx.get_edge_attributes(graph, "weight").items()}
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=labels)

    nxkwargs = {'node_color': nodecolors, 'node_size': 400, 'font_weight': "bold", 
                'font_color': "white", 'cmap': cmap}
    if len(kwargs) > 0:
        for k,v in kwargs.items():
            if k in ["figsize","dpi","facecolor","edgecolor","frameon","clear","layout"]:
                plt.rcParams[f"figure.{k}"] = v
            if k in ["node_color","node_size","font_weight","font_size","font_color","font_weight",
                     "cmap","vmin","vmax","with_labels","arrows","arrowsize","width","ax","hide_ticks"]:
                nxkwargs[k] = v
    nx.draw_networkx(graph, pos, **nxkwargs)

    if len(set(mycolors)) == N and isinstance(mycolors[0], int):
        if np.min(mycolors) == 0 and np.max(mycolors) == N-1:
            all_weights = list(nx.get_edge_attributes(graph, "weight").values())
            if len(all_weights) == 0:
                all_weights = [2.] * N
            headwidth = np.max(np.abs(all_weights)) * 0.02
            for i,j in zip(mycolors[:-1],mycolors[1:]):
                mydist = ( np.array(mypos[j]) - np.array(mypos[i]) ) 
                factor = 0.92
                while np.linalg.norm(mydist*factor) < 30*headwidth and factor > 0.7:
                    factor -= 0.04 
                plt.arrow(mypos[i][0], mypos[i][1], mydist[0]*factor, mydist[1]*factor, 
                      color='r', head_width=headwidth)
            #text = ["(1st)","(2nd)","(3rd)"]
            #text.extend([f"({i+1}th)" for i in range(3, N)])
            #for i in range(N):
            #    plt.annotate(text[i], (pos[i][0]+0.05, pos[i][1]+0.05), size=9, color=nodecolor[i])

    # save / show the graph
    plt.axis("off")
    if filename is not None:
        plt.savefig(filename)
    plt.show()

#######################################################################
# simple plot function 

def plot_histogram(data,                # single dataset or list or array of datasets 
                   data1=None,          # kludge for dealing with 'plot_histogram(xvalues,yvalues,...) 
                   label=None,          # label for each dataset 
                   xvals=None,          # list or array of x-values (same length as datasets)
                   xlabel=None,         # label for x-axis
                   ylabel=None,         # label for y-axis
                   legend_loc=None,     # location of legend box: "best","upper right","upper left", etc)
                   filename=None,       # filename to save the figure
                   **kwargs):           # additional plot parameters (e.g. title=.., figsize=(width,height), etc)

    myarglist = [data1, label, xvals, xlabel, ylabel, legend_loc, filename]
    myargnames = ["label", "xvals", "xlabel", "ylabel", "legend_loc", "filename"]
    myargs = {}
    for n,ptr in zip(myargnames, myarglist[1:]):
        myargs[n] = ptr
    xdata = None
    ydata = data
    if data1 is not None:
        if isinstance(data1, str) or ( isinstance(data1, (list,tuple)) and isinstance(data1[0], str) ):
            for n,ptr in zip(myargnames, myarglist[:-1]):
                if ptr is None:
                    break
                myargs[n] = ptr
        else:
            ydata = data1
            xdata = data        

    if len(kwargs) > 0:
        for k in myargnames:
             if k in kwargs.keys():
                 #print("kwargs:",k,kwargs[k])
                 myargs[k] = kwargs[k]
                 del kwargs[k]
    if xdata is None:
        xdata = myargs["xvals"]
    legend_loc = myargs["legend_loc"] or "best"

    if isinstance(ydata, np.ndarray):
        yvals = ydata.tolist()
    elif isinstance(ydata, dict):
        if xdata is None:
            xdata = list(ydata.keys())
        yvals = list(ydata.values())
    elif isinstance(ydata, (list,tuple)):
        yvals = ydata.copy()
    else:
        print("*** Provide input data in form of list or array or dict")
        return 
    # single or multiple plots?
    if isinstance(yvals[0], (list,tuple)):
        ndim = len(yvals)
        lenarr = []
        for i in range(ndim):
            lenarr.append( len(yvals[i]) )
        nvals = np.max(np.array(lenarr))
    else:
        ndim = 1
        nvals = len(yvals)
        yvals = [yvals]

    myxvals = [None] * ndim
    if xdata is not None: 
        if isinstance(xdata, int):
            # xdata as int (potentially rebinning of data) or list of form (xmin,xmax):
            if xdata == nvals:
                myxvals = list(range(nvals))
            elif nvals%xdata == 0:
                yvals = list(np.array(yvals).reshape((ndim,-1,xdata,nvals//xdata)).mean(-1).mean(1))
                myxvals = list(range(0, nvals, nvals//xdata))
        elif isinstance(xdata, (list,tuple)) and len(xdata)==2:
            myxvals = list(np.linspace(xdata[0], xdata[1], nvals))   
        elif isinstance(xdata, np.ndarray):
            myxvals = xdata.tolist()
        else:
            myxvals = xdata.copy()
        if not isinstance(myxvals[0], (list,tuple)):
            myxvals = [myxvals]

    label = myargs["label"]
    xlabel = myargs["xlabel"]
    ylabel = myargs["ylabel"]
    mytitle = None
    mylabel = [None] * ndim
    if label is not None:
        if isinstance(label, str):
            if ndim == 1:
                # special cases?
                if label.lower().startswith("cost"):
                    mytitle = "Cost function value over iterations"
                    xlabel = xlabel or "Iterations"
                else:
                    mylabel = [label]
            else:
                mytitle = label
        elif isinstance(label, (list,tuple)):
            mylabel = label.copy()
            if len(mylabel) < ndim:
                print("plot_histogram info: not enough labels for all datasets")
                mylabel.extend([None]*(ndim-len(label)))

    # some plt.subplots() settings
    plotargs=[{}] * ndim
    if len(kwargs) > 0:
        if "title" in kwargs.keys():
            mytitle = kwargs["title"]
            del kwargs["title"]
        # kwargs for plt.plot
        plotkwargs = ['data','alpha','color','c','drawstyle','fillstyle',
                      'linestyle','ls','linewidth','lw','marker','markerfacecolor',
                      'markersize','scalex','scaley']
        for k in kwargs.keys():
            if k in plotkwargs:
                val = kwargs[k]
                if not isinstance(val, (list,tuple)):
                    val = [val] * ndim
                for i in range(ndim):
                    plotargs[i][k]=val[i]
        if len(plotargs[0]) > 0:
            for k in plotargs[0].keys():
                del kwargs[k]
        # fmt: '[marker][line][color]' is optional argument of plt.plot
        if "fmt" in kwargs.keys():
            fmt = kwargs["fmt"]
            if isinstance(fmt, str):
                fmt = [fmt] * ndim
            for i,f in enumerate(fmt):
                ls = ""
                for l in list(f):
                    if l in list("bgrcmykw"):
                        plotargs[i]["color"] = l
                    elif l in list(",^ov<>12348spP*hH+xXdD|_"):
                        plotargs[i]["marker"] = l
                    else:
                        ls += l
                if len(ls) > 0:
                    plotargs[i]["linestyle"] = ls
            del kwargs["fmt"]
        if "ax" not in kwargs.keys():
            fig,ax = plt.subplots(**kwargs)
        else:
            ax = kwargs["ax"]
            if isinstance(ax, (list,tuple)):
                ax = np.array(ax)
    else:
        fig,ax = plt.subplots()

    # plot the data
    # multiple subplots? -> each dataset in different subplot
    # or all datasets in a single plot?
    if isinstance(ax, np.ndarray):
        axes = ax.flatten().tolist()
        if myxvals[0] is None:
            for i in range(ndim):
                myxvals[i] = list(range(0, lenarr[i]))
        else:
           if len(myxvals) < ndim:
               for j in range(len(myxvals), ndim):
                   myxvals.append(myxvals[0])
        if xlabel is None:
            xlabel = [None] * ndim
        elif not isinstance(xlabel, (list,tuple)):
            if ax.ndim > 1:
                xlabel = [None]*ax.shape[0]*(ax.shape[1]-1) + [xlabel]*ax.shape[0]
            else:
                xlabel = [xlabel]*ax.shape[0]

        for j,(x,y,iax) in enumerate(zip(myxvals, yvals, axes)):
            iax.plot(x, y, label=mylabel[j], **plotargs[j])
            iax.set_xlabel(xlabel[j])
            iax.set_ylabel(ylabel)
            if mylabel[j] is not None:
                iax.legend(fontsize="medium", loc=legend_loc)
        if mytitle is not None:
            fig = axes[0].get_figure()
            fig.suptitle(mytitle, fontweight="bold", fontsize="large")
    else:
    # single plot with all datasets
        if ndim > 1 and not np.all(np.array(lenarr)==nvals):
            imax = int(np.argmax(np.array(lenarr)))
            if myxvals[0] is not None:
                if len(myxvals[0]) < nvals:
                    myxvals[0] = myxvals[imax]
            for i in range(ndim):
                if lenarr[i] < nvals:
                    #yvals[i].extend( [yvals[i][-1]]*(nvals-lenarr[i]) )
                    yvals[i].extend( [None]*(nvals-lenarr[i]) )
        
        for i in range(ndim):
            if myxvals[i] is None:
                ax.plot(yvals[i], label=mylabel[i], **plotargs[i])
            else:
                ax.plot(myxvals[i], yvals[i], label=mylabel[i], **plotargs[i])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if mylabel[0] is not None:
            plt.legend(fontsize="medium", loc=legend_loc)
        if mytitle is not None:
            plt.title(mytitle, fontweight="bold", fontsize="large")
        
    if isinstance(ax, np.ndarray) and xlabel[0] is None:
        plt.tight_layout() 
    if myargs["filename"] is not None:
        plt.savefig(myargs["filename"])
    plt.show()

#######################################################################

def plot_probabilities(probabilities,    # dict with key=state or list of probabilities for all possible states
                       all_states=False, # all possible states on x-axis if True
                       filename=None,    # save figure in file
                       **kwargs):        # additional plot parameters
    """
  plot_probabilities( probabilities, all_states, filename, **kwargs )
    create histogram with probability (or counts) entries
    input:  probabilities: dict with state as key and probability (or counts) as value or 
                           list (or array) of probabilities for all possible states   (dict or array or list) 
            all_states:    =True: plot for all possible states; =False (default): plot only states with prob.>0 
            filename:      save figure in file                              (string, optional)
            kwargs:        additional plot parameters                 (dict entries, optional)
    """
    mydict = {}
    keylen = 1
    if isinstance(probabilities, dict):
        for k,v in probabilities.items():
            if isinstance(k, int):
                k = bin(k)[2:]
            if len(k) > keylen:
                keylen = len(k)
            if v > 0.0:
                mydict[k] = v
    else:
        # list or array of samples?
        samples = []
        tmp = probabilities.tolist() if isinstance(probabilities, np.ndarray) else probabilities
        if isinstance(tmp[0], (list,tuple)):    # list of values per sample, e.g. [0,1,1,0,1]
            for vals in tmp:
                samples.append("".join([ str(i) for i in vals ]))
        elif isinstance(tmp[0], str):           # assume sample like "01101"
            for vals in tmp:
                samples.append( vals )
        elif isinstance(tmp[0], int):
            for vals in tmp:
                samples.append( bin(vals)[2:] )

        if len(samples) > 0:
            mydict = defaultdict(int)
            for entry in samples:
                 mydict[entry] += 1
                 if len(entry) > keylen:
                     keylen = len(entry)
        else:
        
            if len(tmp)>0: 
                keylen = math.ceil(math.log2(len(tmp)))
            for j in range(len(tmp)):
                if probabilities[j] > 0.0:
                    mydict[ bin(j)[2:] ] = probabilities[j]

    myfigsize=(8,4)
    if all_states:
        xvalues = [ bin(j)[2:].zfill(keylen) for j in range(2**keylen) ]
        probs = [ 0.0 for _ in range(2**keylen) ]
        for k,v in mydict.items():
            probs[ int(k,2) ] = v
        if keylen > 5 and ( len(kwargs) == 0 or
             ("figsize" not in kwargs.keys() and "ax" not in kwargs.keys()) ):
            myfigsize=(16,8)
    else:
        # sort by key
        xvalues = []
        probs = []
        for k in sorted(mydict):
            xvalues.append( k.zfill(keylen) )
            probs.append( mydict[k] )

    tmpdict = {}
    for arg in ["title", "xlabel", "ylabel", "label", "legend_loc"]:
        tmpdict[arg] = None
        if len(kwargs) > 0:
            if arg in kwargs.keys():
                tmpdict[arg] = kwargs[arg]
                del kwargs[arg]
    if len(kwargs) > 0: 
        if "ax" in kwargs.keys():
            ax = kwargs("ax")
        else:
            if "figsize" not in kwargs.keys():
                kwargs["figsize"] = myfigsize
            fig, ax = plt.subplots(**kwargs)
    else:
        fig,ax = plt.subplots(figsize=myfigsize)

    # plot the data
    ax.bar(xvalues, probs, label=tmpdict["label"])
    ax.set_xticks(range(len(xvalues)), xvalues)
    ax.set_xticklabels(xvalues, rotation=90)
    ax.set_xlabel(tmpdict["xlabel"])
    ax.set_ylabel(tmpdict["ylabel"])
    ax.margins(x=0)
    plt.title(tmpdict["title"], fontweight="bold", fontsize="large")
    if tmpdict["label"] is not None:
       loc = tmpdict["legend_loc"] or "best"
       plt.legend(tmpdict["label"], loc=loc)
    if filename is not None:
        plt.savefig(filename)
    #if "ax" not in kwargs.keys():
    plt.show()

##############################################################################
#### helper functions for Traveling Salesman Problem
#
# generate coefficients for TSP starting at location 0 using (N-1)**2 qubits 
#      (in matrix form: rows=locations, columns=time steps; ideally initialized with superposition of W-states)
#  (based on "Important Quantum Gates for Quantum Algorithms of Travelling Salesman Problem",
#    by T.J.H.Sinaga et al, in Proceedings of ICoABCD 2023)
#
# tsp_generate_coeffs(J_or_G,           # Ising matrix (distance matrix) or Graph 
#                     scale=2*np.pi,    # scale factor if scale<0 or normalize to scale value
#                     epsilon=0.0,      # set coefficient to zero if scaled coefficient < epsilon 
#                     **kwargs)         # other parameters (not used yet)
#
##################################################################

import numpy as np
import sympy as sp
import networkx as nx
import math

def tsp_generate_coeffs(J_or_G, scale=2*np.pi, epsilon=0.0, **kwargs):
# input J_or_G: Ising matrix (distance matrix) or Graph

    if isinstance(J_or_G, nx.classes.graph.Graph):
        J = nx.to_numpy_array(J_or_G)
    elif isinstance(J_or_G, list) and isinstance(J_or_G[0], list):
        J = np.array(J_or_G)
    elif isinstance(J_or_G, np.ndarray):
        J = J_or_G
    else:
        print("Provide Ising matrix (distance matrix) or Graph object")
        return None

    N = J.shape[0]
    # matrix for all possible edges to implement constraints
    X = np.empty((N,N), dtype=object)
    for ii in range(1, N):
        for jj in range(1, N):
            X[ii,jj] = sp.Symbol(f"x_{ii}{jj}", bool=True)
    # starting at node 0:
    X[0,:] = 0
    X[:,0] = 0
    X[0,0] = 1
    # reduced vector for PauliZ operations
    # scratch vector to hold all variables of matrix X
    # with V_i = 1/2 (I - Z_i) 
    Z = np.empty(((N-1)**2), dtype=object)
    for ii in range(len(Z)):
        Z[ii] = sp.Symbol(f"z_{ii}", bool=True)
    V = np.empty_like(Z)
    for ii in range(len(V)):
        V[ii] = sp.Symbol(f"v_{ii}", bool=True)

    # multiplier (weight for constraint terms)
    A = 10 * np.max(J)
    B = 1             # edge weights (distance)

    # Hamiltonian 
    # constraint term1: one node at a time
    term1 = sp.Integer(0)
    for i in range(N):
        term1 += sp.expand((1 - np.sum(X[i,:])) ** 2)

    # constraint term2: visit each node once
    term2 = sp.Integer(0)
    for j in range(N):
        term2 += sp.expand((1 - np.sum(X[:,j])) ** 2)

    # constraint term3: distance between nodes
    term3 = sp.Integer(0)
    for u in range(N):
        for v in range(u+1, N):
            term3_partial = 0
            for i in range(N-1):
                term3_partial += X[u,i] * X[v,i+1]
            term3_partial += X[u,N-1] * X[v,0]
            term3 += float(J[u,v]) * term3_partial

    H = sp.expand(A*term1 + A*term2 + B*term3)

    # replace squares as V_i is projection operator (V_i=V_i*V_i)
    for xr in X[1:,1:]:
        for x in xr:
            H = H.subs(x**2, x)

    # transform matrix to reduced vector
    for i, xr in enumerate(X[1:,1:]):
        for j,x in enumerate(xr):
            H = H.subs(x, V[i+j*(N-1)])

    # replace V_i=1/2 (I-Z_i)
    for i, v in enumerate(V):
        H = H.subs(v, 1 / 2 * (1 - Z[i]))
    H = H.expand()

    Hterms = sp.Add.make_args(H)
    coeff_terms = {}
    for v in Hterms:
        vv = str(v).split("*")
        # check whether first string is a number
        tmp = vv[0].replace(".","",1).replace("e-","",1).replace("e","",1)
        if not tmp.lstrip("-").isdigit():
            vv.insert(0, "-1" if tmp.startswith("-") else "1")
        if len(vv) == 1:
            key = "offset"
        else:
            key = "*".join(vv[1:])
        coeff_terms.update({key: float(vv[0])})

    # scale (discard "offset" for scaling)
    maxval = np.max(np.abs(list(coeff_terms.values())[1:]))
    scale = scale/maxval if scale > 0. else abs(scale)
    for k,v in coeff_terms.items():
        if abs(v*scale) > epsilon:
            coeff_terms[k] = round(v*scale,3)

    return coeff_terms

###########

def tsp_get_most_probable_paths(probabilities, nresults=5, tsp_fidelity=2, list_states=False, rounding_digits=4):
    """
 tsp_get_most_probable_paths(probabilites,      # dict with key=state or list of probabilities for all possible states
                             nresults=3,        # number of N best results to report in output
                             tsp_fidelity=2,    # output tsp_matrix
                             list_states=False, # list of all states contributing to this output (see below)
                             rounding_digits=4) # rounding of probabilities for output
     returns sorted dict of results:
       for tsp_fidelity=0: all results independent of whether they represent a 'good' path
       for tsp_fidelity=1: results where each place was visited (even if 2 places are visited at the same time)
       for tsp_fidelity=2: 'good' results (which contain all places (nodes) visited at different times).
     output:
       if list_states==False (default):
          dict of form {'good_path': [colorlist, probability]}
       if list_states==True:
          dict of form {'good_path': [colorlist, probability, list_of_states]}
        (note: good path is e.g. ['1000','0001','0100','0010'] but also ['1000','0101','0110','0001'] (where
              first all unique timestamps are counted, the remaining are split between remaining
              places (nodes) if they have a '1' for these timestamps: (here: ['1000','0100','0010','0001']))
    """
    maxstates = len(probabilities)
    arr = get_most_probable_states(probabilities, maxstates, rounding_digits=rounding_digits)

    myarr = {}
    mylen = len(list(arr.keys())[0])
    Nvals = mylen//math.ceil(math.sqrt(mylen))  # number of places (nodes) without starting point 

    for k,v in arr.items():
        tmp = np.array([int(i) for i in k]).reshape((Nvals,Nvals))

        ### is checking for rank a good idea?
        ### e.g. [[0,1,1,0],[1,0,0,0],[0,1,1,0],[0,0,0,1]] can be regarded as valid
        ###  but [[0,0,1,1],[0,0,0,0],[1,0,0,0],[0,1,1,0]] not
        #if np.linalg.matrix_rank(tmp) < Nvals:
        #    continue
        # np.any(array&condition, axis=1): returns True/False for each row
        if tsp_fidelity > 0 and np.sum(np.any(tmp == 1, axis=1)) < Nvals:
             continue
        mycolors = [ 0 for _ in range(Nvals) ]
        tt = np.where(tmp == 1)              # tt[0]=array with '1' in rows, tt[1]=array with '1' in columns
        uu = np.array(np.unique(tt[0], return_counts=True))   # count occurrence of row number in tt[0]
        dd = []                                               # save row numbers with multiple '1'
        for i,j in zip(uu[0],uu[1]):         # uu[0]=row number, uu[1] number of '1'
            if j == 1:
                ll = np.where(tt[0]==i)[0][0]
                mycolors[i] = tt[1][ll]+1
            else:
                dd.append([i, np.argwhere(tt[0]==i).T.tolist()[0] ])
        for i,p in dd:
            m = [ tt[1][j]+1 for j in p ]
            nn = [ n for n in m if n not in mycolors ]
            if len(nn) > 0:
                mycolors[i] = nn[0]
        # skip all cases where no valid order is found (each number once)
        mycolors = [0] + mycolors
        if tsp_fidelity == 2 and sorted(mycolors) != list(range(Nvals+1)):
            continue
        mystate = [ "".join([str(i) for i in irow]) for irow in tmp.tolist() ]
        myarr[k] = [v, mystate, mycolors]

    colsdict = defaultdict(int)
    if list_states:
        statedict = defaultdict(list)
    for k,val in myarr.items():
        mykey = "".join([ str(i) for i in val[2] ])
        colsdict[mykey] += val[0]
        if list_states:
            statedict[mykey].append(val[1])
    sorted_cols = OrderedDict(sorted(colsdict.items(), key=lambda item: item[1], reverse=True))
    mydict = {}
    for k in list(sorted_cols)[:nresults] :
        if list_states:
            mydict[k] = [ [int(i) for i in list(k)], round(sorted_cols[k], rounding_digits), statedict[k] ]
        else:
            mydict[k] = [ [int(i) for i in list(k)], round(sorted_cols[k], rounding_digits) ]
    return mydict
    
