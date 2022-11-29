function [Fe] = COMPUTE_Fe_FORCE(f,node, elementCoor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
eps = [-1/sqrt(3) 1/sqrt(3)];
ws = [1 1];
numberNodes = size(node,1);
Fe = zeros ( numberNodes , 1 );
for g=1:2
    Ne = 1/2 * [(1-eps(g)) (1+eps(g))];
    xe = Ne * elementCoor';
    q = Ne' * subs(f,xe);
    Fe = Fe + ws(g) * q;
end
end

