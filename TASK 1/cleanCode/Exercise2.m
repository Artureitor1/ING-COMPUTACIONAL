%APARTADO 2
clc
clear all 
clc
syms theta u(x) 

L = 1;
g = 0.01;
rho = pi^2/L^2;
s = g*rho^2;
F_AE = g*pi^2/L;



du = diff(u,x);
ode = diff(u,x,2)+u*rho-s*x^2 == 0;
usol(x) = dsolve(ode, u(0) == -g, du(L) == F_AE); 
disp(simplify(usol));
%PLot
figure
hold on
title('Graph of displacement along the bar by exact method', 'Interpreter','latex', 'FontSize',16)
fplot(usol,[0 L],'Linewidth',2)
xlabel('x', 'Interpreter','latex', 'FontSize',15)
ylabel('displacement', 'Interpreter','latex', 'FontSize',15)
grid on
grid minor
hold off
