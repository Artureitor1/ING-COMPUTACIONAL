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
elementsNumber = 40;
[nodes, coords] = ShapeFunctionsFiniteElement1D(elementsNumber, 0, L);
nodes = 1:length (coords);
restringedNodes = [1 -g];
restringedForce= [nodes(end) F_AE];
conectivityMatrix = computeConectivityMatrix1D(nodes);
%Solver 
nodeDisplacement = computeDisplacement(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho); 

%Exact Solution:
syms x
uExact = (61740346358469708232516810123*cos((2778046668940015^(1/2)*x)/16777216))/6174034635846970641698934560180 + (8773830921664641*x^2)/88897493406080480 + (2778046668940015^(1/2)*sin((2778046668940015^(1/2)*x)/16777216)*(66293192113335022366330415812131684352*2778046668940015^(1/2)*sin(2778046668940015^(1/2)/16777216) - 10977124066582370799746477992578602748793214053))/(18416558152385082682136054268153845979438084432645324800*cos(2778046668940015^(1/2)/16777216)) - 154350865896174268311882694656/7717543294808713302123668200225;

%Plot 
figure
hold on
title('Graph of displacement along the bar by exact method and FEM method', 'Interpreter','latex', 'FontSize',16)
fplot(uExact,[0 L],'Linewidth',2)
plot(coords,nodeDisplacement,'Linewidth',2);
grid on
grid minor
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('displacement', 'Interpreter','latex', 'FontSize',15)
legend('Exact solution', 'FEM solution', 'location', 'northeast', 'Interpreter', 'latex', 'FontSize',12)
hold off


