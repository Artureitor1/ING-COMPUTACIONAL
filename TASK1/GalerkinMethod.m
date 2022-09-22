function [u]= GalerkinMethod(N,f,b,g,rho) ;
syms x % 
N_1 = subs(N,1) ;
B = diff(N,x); 
NtN = N.'*N;
BtB = B.'*B;  
K = int(BtB-rho*NtN,0,1); 
F= int(N.'*f,0,1) + N_1.'*b; 
r = 1; 
l = 2:length(N); 
dl = K(l,l)\(F(l)-K(l,r)*g); 
u  = g + N(l)*dl; 
end