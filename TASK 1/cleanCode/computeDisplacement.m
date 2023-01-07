function [d] = computeDisplacement(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho); 

stiffnesMatrix = AssemblyK(coords, conectivityMatrix, rho);

Force = AssemblyF(coords, conectivityMatrix, f, nodes);
ft = zeros(size(Force));
ft(restringedForce(:,1)) = restringedForce(:,2);
Force = Force + ft;
freeNodes = setdiff(nodes,restringedNodes(:,1));
dr = zeros(length(nodes),1);
dl = zeros(length(nodes),1);
dr(restringedNodes(:,1)) = restringedNodes(:,2);
dl(freeNodes) = stiffnesMatrix(freeNodes,freeNodes) \ (Force(freeNodes)-stiffnesMatrix(freeNodes,restringedNodes(1,1))*restringedNodes(1,2));
d = dl+dr;
