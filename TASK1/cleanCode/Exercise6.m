%Apartado 6
clc
clear all
close all


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
elementsNumber = 30;
[nodes, coords] = ShapeFunctionsFiniteElement1D(elementsNumber, 0, L);
nodes = 1:length (coords);
restringedNodes = [1 -g];
restringedForce= [nodes(end) F_AE];
freeNodes = setdiff(nodes,restringedNodes(1,1));
conectivityMatrix = computeConectivityMatrix1D(nodes);
%Solver 
[uGalerking, d] = GalerkinMethod(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho); 

hold on
xPlot = (0:0.01:L);
plot(coords,d)
grid on 
hold off


