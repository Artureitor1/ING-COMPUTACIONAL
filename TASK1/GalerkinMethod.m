function [u, d]= GalerkinMethod(N,f,b,g,rho, CN, coords)
syms x % 


K = AssemblyK(coords, CN, rho );
F= AssemblyF(coords, CN, f, N);
ft = zeros(size(F));
ft(end) = b;
F = F + ft;
r = 1; 
l = 2:length(N); 
dl = K(l,l) \ (F(l)-K(l,r)*g);
d = [g;dl];
u  = +g + N(l)*dl; 
end