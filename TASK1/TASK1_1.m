%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear a
%APARTADO 1


syms theta u(x) 

L = 1;
g = 0.01;
rho = pi^2/L^2;s = g*rho^2;
f =  -s*x^2;
nelms = 10;
[N, coords] = ShapeFunctionsFiniteElement1D(nelms, 0, L);
CN = [uint32(1):uint32(size(coords,2)-1); uint32(2):uint32(size(coords,2))]';
figure(2)
du = diff(u,x);
usol(x) = dsolve(diff(u,2) + u*rho == -f, u(0) == -g, du(L) == g*pi^2/L); 
%[uaprox] = MinimizationResidual(N,f,g*pi^2/L,-g);
hold on
fplot(usol,[0 L]);


%APARTADO 2

hold on
[uGalerking, d] = GalerkinMethod(N,f,g*pi^2/L,-g,rho, CN, coords); 
%fplot(uGalerking,[0,1]);
plot(coords, d)
hold off


