function [u, d] = GalerkinMethod(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho); 

syms x 
stiffnesMatrix = AssemblyK(coords, conectivityMatrix, rho);

Force = AssemblyF(coords, conectivityMatrix, f, nodes);
ft = zeros(size(Force));
ft(restringedForce(1,1)) = restringedForce(1,2);
Force = Force + ft;

r = 1; 
l = 2:length(nodes); 
dl = stiffnesMatrix(l,l) \ (Force(l)-stiffnesMatrix(l,r)*restringedNodes(1,2));
d = [restringedNodes(1,2);dl];
u  = +restringedNodes(1,2) + nodes(l)*dl; 
