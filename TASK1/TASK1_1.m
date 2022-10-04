%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%APARTADO 1
syms theta u(x) 

L = 1;
g = 0.01;
numElem = 5;
rho = pi^2/L^2;
s = g*rho^2;
f =  -s*x^2;
N = [1 x x^2 x^3 x^4];

du = diff(u,x);
usol(x) = dsolve(diff(u,2) + u*rho == -f, u(0) == -g, du(L) == g*pi^2/L); 
hold on
fplot(usol,[0 L]);

uGalerking = GalerkinMethod(N,f,g*pi^2/L,-g,rho); 
fplot(uGalerking,[0,1]);
hold off

%APARTADO 2

COOR = linspace(0,L,numElem);





