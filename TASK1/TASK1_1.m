%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%APARTADO 1

syms theta u(x) 

L = 1;
g = 0.01;
rho = pi^2/L^2;
s = g*rho^2;



du = diff(u,x);
ode = diff(u,x,2)+u*rho-s*x^2 == 0;
usol(x) = dsolve(ode, u(0) == -g, du(L) == g*pi^2/L); 
fplot(usol,[0 L]);
