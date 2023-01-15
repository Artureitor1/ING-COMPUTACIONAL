import numpy as np
import matplotlib.pyplot as plt

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

def computeDisplacement(nodes, restringedNodes, restringedForce, coords, conectivityMatrix, f, F_AE, rho):

    stiffnesMatrix = AssemblyK(coords, conectivityMatrix, rho);

    Force = AssemblyF(coords, conectivityMatrix, f, nodes);
    ft = np.zeros(len(Force));
    ft[int(restringedForce[0])] = restringedForce[1]
    Force = Force + ft
    dr = np.zeros(len(nodes))
    dl = np.zeros(len(nodes))
    dr[int(restringedNodes[0])] = restringedNodes[1]
    sub_stiffnesMatrix = stiffnesMatrix[1:,1:]
    sub_stiffnesMatrix_2 = stiffnesMatrix[1:,int(restringedNodes[0])]
    sub_Force = Force[1:]
    dl[1:] = np.dot(np.linalg.inv(sub_stiffnesMatrix),(sub_Force-sub_stiffnesMatrix_2*restringedNodes[1]))
    d = dl+dr;
    return d

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

# Problem data
L = 1.0
g = 0.01
rho = (np.pi / L)** 2
s = g * rho**2
f = lambda x: -s * x**2
F_AE = (g * np.pi**2) / L

elementsNumber = 40

coords = ShapeFunctionsFiniteElement1D(elementsNumber, 0, L);
nodes = range(len(coords))

restringedNodes = np.array([0, -g])
restringedForce = np.array([nodes[-1], F_AE])

conectivityMatrix = computeConectivityMatrix1D(nodes)

nodeDisplacement = computeDisplacement(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho); 

plt.plot(coords, nodeDisplacement, 'o-')
plt.show()