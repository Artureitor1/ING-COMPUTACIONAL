%Apartado 4
clc
clear all
close all

syms x

L = 1;
g = 0.01;
rho = (pi/L)^2;
s = g*rho^2;
f =  -s*x^2;
F_AE = (g*pi^2)/L;
N = [1 x x^2 x^3 x^4];
B = diff(N,x);
d = zeros(length(N),1);
nodes = (1:length(N));
restringedNodes = 1;
freeNodes = setdiff(nodes,restringedNodes);

%WeakForm Solver:
stiffnessMatrix = int(B.'*B-rho*(N.'*N),x,[0 L]);
force = int(N.'*f,x,[0,L])+subs(N.',L)*F_AE;

d(restringedNodes) = -g;
d(freeNodes) = stiffnessMatrix(freeNodes,freeNodes)\(force(freeNodes)-stiffnessMatrix(freeNodes,restringedNodes)*d(restringedNodes));

%Plot 
hold on
xPlot = (0:0.01:L);
plot(xPlot,polyval(flipud(d),xPlot));
grid on 
hold off




