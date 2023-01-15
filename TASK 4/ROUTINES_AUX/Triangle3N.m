function  [NeXi, BeXi] = Triangle3N(xi)
% Shape functions and derivatives for 2-node linear element 
 
% Matrix of shape functions
e1 = xi(1);     e2 = xi(2);     

NeXi = [1-e1-e2 e1 e2];
% Matrix of the gradient of shape functions 
BeXi = [1 0 0;
        0 1 0;
        0 0 1];