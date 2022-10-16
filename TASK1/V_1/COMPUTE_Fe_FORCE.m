function [outputArg1] = COMPUTE_Fe_FORCE(f,N, COOR_e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
eps = [-1/sqrt(3) 1/sqrt(3)];
ws = [1 1];
nnode = size(N,1);
Fe = zeros ( nnode , 1 );
for g=1:2
    Ne = 1/2 * [(1-eps(g)) (1+eps(g))];
    xe = Ne * COOR_e';
    q = Ne' * subs(f,xe);
    Fe = Fe + ws(g) * q;
end

outputArg1 = Fe;
end

