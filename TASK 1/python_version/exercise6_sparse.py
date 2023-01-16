import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sparse
import time
from random import random


def COMPUTE_Fe_FORCE(f,node, elementCoor):
    eps = [-1/np.sqrt(3), 1/np.sqrt(3)];
    ws = [1, 1];
    Fe = 0.0
    for g in range(len(ws)):
        Ne = 0.5 * np.array([1-eps[g], 1+eps[g]]);
        xe = np.dot(Ne,np.array(elementCoor).T);
        q = Ne * f(xe);
        Fe = Fe + np.dot(ws[g],q);

    return Fe

def AssemblyK(coor, conectivityMatrix, rho):
    numberElements = len(conectivityMatrix) 
    numberNodes = int(conectivityMatrix[-1][-1]+1);
    nodePerElement = len(conectivityMatrix[0]) ;
    stiffnesMatrix =np.zeros((numberNodes, numberNodes));
    for element in range(numberElements): 
        elementNodes = conectivityMatrix[element]    # Global numbering of nodes of element "e"
        elementCoor = [coor[i] for i in elementNodes]
        he = elementCoor[1]-elementCoor[0] # Size finite element
        # Elemental matrix
        Ke = 1/he* np.array([[1,-1],[-1,1]]) - rho * he/2 * np.array([[2/3, 1/3], [1/3, 2/3]])
        # Assembly
        for a in range(nodePerElement):
            for b in range(nodePerElement):
                A = conectivityMatrix[element][a];
                B = conectivityMatrix[element][b];
                stiffnesMatrix[A][B] = stiffnesMatrix[A][B]+ Ke[a][b];
    
    return stiffnesMatrix  

def AssemblyF(coor, conectivityMatrix, f, nodes):
    numberElements = len(conectivityMatrix) 
    numberNodes = int(numberElements+1);
    nodePerElement = len(conectivityMatrix[0]) ;
    Force = np.zeros(numberNodes);
    for element in range(numberElements): 
        elementNodes = conectivityMatrix[element]    # Global numbering of nodes of element "e"
        elementCoor = [coor[i] for i in elementNodes]
        he = elementCoor[1]-elementCoor[0] # Size finite element
        # Elemental vector
        fe = np.dot((he/2),COMPUTE_Fe_FORCE(f, nodes, elementCoor))
        # Assembly
        for a in range(nodePerElement):
            A = conectivityMatrix[element,a];
            Force[A] = Force[A]+ fe[a];

    return Force

def computeDisplacement(nodes, restringedNodes, restringedForce, coords, conectivityMatrix, f, F_AE, rho, time_stamps):

    time_stamps['assembly_start'] = float(time.perf_counter() - time_stamps['start'])
    stiffnesMatrix = sparse.csr_matrix(AssemblyK(coords, conectivityMatrix, rho))
    time_stamps['assembly_end'] = float(time.perf_counter() - time_stamps['start'])

    time_stamps['force_start'] = float(time.perf_counter() - time_stamps['start'])
    Force = sparse.csr_array(AssemblyF(coords, conectivityMatrix, f, nodes))
    time_stamps['force_end'] = float(time.perf_counter() - time_stamps['start'])

    ft = np.zeros((Force.shape[1]));
    ft[int(restringedForce[0])] = restringedForce[1]
    Force = Force + ft
    dr = np.zeros(len(nodes))
    dl = np.zeros(len(nodes))
    dr[int(restringedNodes[0])] = restringedNodes[1]
    sub_stiffnesMatrix = stiffnesMatrix[1:,1:]
    sub_stiffnesMatrix_2 = stiffnesMatrix[1:,int(restringedNodes[0])]
    sub_Force = Force[:,1:].T
    time_stamps['solve_start'] = float(time.perf_counter() - time_stamps['start'])
    # = 
    dl[1:] = sparse.linalg.spsolve(sub_stiffnesMatrix,(sub_Force-sub_stiffnesMatrix_2*restringedNodes[1]))
    d = dl+dr;
    time_stamps['solve_end'] = float(time.perf_counter() - time_stamps['start'])
    return d, time_stamps

def ShapeFunctionsFiniteElement1D(n, x0, x1):
    nnodes = n + 1
    coords = np.linspace(x0, x1, nnodes)

    return coords

def computeConectivityMatrix1D(nodes):
    conectivityMatrix = np.zeros((len(nodes)-1,2), dtype=int);
    for i in range(len(nodes)-1):
        conectivityMatrix[i,0] = int(i)
        conectivityMatrix[i,1] = int(i+1)
    return conectivityMatrix

def computeDisplacementAll(L, g, elementsNumber, time_stamps):

    rho = np.power((np.pi**1)/L, 2)
    s = g * np.power(rho, 2)
    f = lambda x: -s * np.power(x,2)
    F_AE = (g * np.pi**2) / L
    time_stamps['shape_func_start'] = float(time.perf_counter() - time_stamps['start'])
    coords = ShapeFunctionsFiniteElement1D(elementsNumber, 0, L);
    time_stamps['shape_func_end'] = float(time.perf_counter() - time_stamps['start'])
    nodes = range(len(coords))

    restringedNodes = np.array([0, -g])
    restringedForce = np.array([nodes[-1], F_AE])

    time_stamps['conn_start'] = float(time.perf_counter() - time_stamps['start'])
    conectivityMatrix = computeConectivityMatrix1D(nodes)
    time_stamps['conn_end'] = float(time.perf_counter() - time_stamps['start'])

    time_stamps['disp_start'] = float(time.perf_counter() - time_stamps['start'])
    nodeDisplacement, time_stamps = computeDisplacement(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho, time_stamps); 
    time_stamps['disp_end'] = float(time.perf_counter() - time_stamps['start'])
    return nodeDisplacement, coords, time_stamps

def give_times(n = 15, el = 40):
    times = []
    for _ in range(n):
        _, _, time_stamps = computeDisplacementAll(1.0, 0.01, el, {'start': time.perf_counter()})
        time_stamps['end'] = float(time.perf_counter() - time_stamps['start'])
        time_stamps['start'] = 0
        times.append(time_stamps)
    return times
    
if __name__ == "__main__":
    n_elms = [1, 10, 100, 1000, 2000, 3000]
    results = []
    for el in n_elms[0:2]:
        times = give_times(n=10, el = el)
        total = np.mean([t['end'] - t['start'] for t in times])
        results.append(total)

        print("Total time: ", total)
        print("Solve time: ", np.mean([t['solve_end'] - t['solve_start'] for t in times]))
        print("Assembly time: ", np.mean([t['assembly_end'] - t['assembly_start'] for t in times]))
        print("Force time: ", np.mean([t['force_end'] - t['force_start'] for t in times]))
   
    plt.plot(n_elms, results)
    plt.show()

    
    file = open("results_sparse.txt", "w")
    for i in range(len(n_elms)):
        file.write(str(n_elms[i]) + " " + str(results[i]) + "\n")
    file.close()

