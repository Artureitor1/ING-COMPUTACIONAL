clear all
close all
elementsNumberArray = [4, 6, 8, 10];
errorArray= zeros(size(elementsNumberVector,2),1);

for n = 1:size(nelms,2)
    elementNumber= elementsNumberArray(n);
    [nodeDisplacement,conectivityMatrix,coords]  = computeFEMmethod1D(elementNumber);
    error = computeError(nodeDisplacement); 
    errorArray(n) = error;
end


%Functions
function [nodeDisplacement,conectivityMatrix,coords] = computeFEMmethod1D (elementsNumber)
    %% Data
    %Problem Data
    syms x
    L = 1;
    g = 0.01;
    rho = (pi/L)^2;
    s = g*rho^2;
    f =  -s*x^2;
    F_AE = (g*pi^2)/L;
    
    %Numeric Data
    [nodes, coords] = ShapeFunctionsFiniteElement1D(elementsNumber, 0, L);
    nodes = 1:length (coords);
    restringedNodes = [1 -g];
    restringedForce= [nodes(end) F_AE];
    conectivityMatrix = computeConectivityMatrix1D(nodes);
    %Solver
    nodeDisplacement = computeDisplacement(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho);
end
function error = computeError(nodeDisplacement)
end 