import numpy as np

def COMPUTE_Fe_FORCE(f,node, elementCoor):
    eps = [-1/np.sqrt(3), 1/np.sqrt(3)];
    ws = [1, 1];
    Fe = 0.0
    for g in range(len(ws)):
        Ne = 0.5 * np.array([1-eps[g], 1+eps[g]]);
        xe = Ne * np.array(elementCoor).T;
        q = Ne * f(xe);
        Fe = Fe + ws[g] * q;

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
    numberNodes = int(conectivityMatrix[-1][-1]+1);
    nodePerElement = len(conectivityMatrix[0]) ;
    Force = np.zeros(numberNodes);
    for element in range(numberElements): 
        elementNodes = conectivityMatrix[element]    # Global numbering of nodes of element "e"
        elementCoor = [coor[i] for i in elementNodes]
        he = elementCoor[1]-elementCoor[0] # Size finite element
        # Elemental vector
        fe = (he/2) * COMPUTE_Fe_FORCE(f, nodes, elementCoor)
        # Assembly
        for a in range(nodePerElement):
            A = conectivityMatrix[element,a];
            Force[A] = Force[A]+ fe[a];
    return Force

def computeDisplacement(nodes, restringedNodes, restringedForce, coords, conectivityMatrix, f, F_AE, rho):

    stiffnesMatrix = AssemblyK(coords, conectivityMatrix, rho);

    Force = AssemblyF(coords, conectivityMatrix, f, nodes);
    ft = np.zeros(len(Force));
    ft[restringedForce[:,1]] = restringedForce[:,2]
    Force = Force + ft
    dr = np.zeros(len(nodes))
    dl = np.zeros(len(nodes))
    dr[restringedNodes[:,1]] = restringedNodes[:,2]
    stiffnesMatrix[1:,1:] /(Force[1:]-stiffnesMatrix[1:,restringedNodes[1,1]]*restringedNodes[0,1])
    d = dl+dr;
    return d

def ShapeFunctionsFiniteElement1D(n, x0, x1):
    nnodes = n + 1
    coords = np.linspace(x0, x1, nnodes)
    N = []
    for i in range(nnodes):
        if i == 0:
            x1, x2 = coords[i], coords[i+1]
            N.append(lambda x: (x2-x)/(x2-x1))
        
        elif i == nnodes-1:
            x1, x2 = coords[i-1], coords[i]
            N.append(lambda x: (x-x1)/(x2-x1))
        
        else:
            x1, x2, x3 = coords[i-1], coords[i], coords[i+1]
            N.append(lambda x: (x-x1)/(x2-x1) if x < x2 else (x3-x)/(x3-x2))

    return N, coords

def computeConectivityMatrix1D(nodes):
    conectivityMatrix = np.zeros((len(nodes)-1,2), dtype=int);
    for i in range(len(nodes)-1):
        conectivityMatrix[i,0] = int(i)
        conectivityMatrix[i,1] = int(i+1)
    return conectivityMatrix

# Problem data
L = 1.0
g = 0.01
rho = np.power(np.pi / L, 2)
s = g * np.power(rho, 2)
f = lambda x: -s * np.power(x,2)
F_AE = (g * np.power(np.pi, 2)) / L

elementsNumber = 40
nodes, coords = ShapeFunctionsFiniteElement1D(elementsNumber, 0, L);
nodes = range(1, len(coords))
restringedNodes = [1, -g]
restringedForce = [nodes[-1], F_AE]
conectivityMatrix = computeConectivityMatrix1D(nodes)

nodeDisplacement = computeDisplacement(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho); 
print("ok")