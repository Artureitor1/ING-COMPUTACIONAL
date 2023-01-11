function  [NeXi, BeXi] = Hexahedra8N(xi)
% Shape functions and derivatives for 2-node linear element 
 
% Matrix of shape functions
e1 = xi(1);     e2 = xi(2);     e3 = xi(3);
p1 = 1+e1;      p2 = 1+e2;      p3 = 1+e3;
n1 = 1-e1;      n2 = 1-e2;      n3 = 1-e3;

NeXi =(1/8)*[n1*n2*n3 p1*n2*n3 p1*p2*n3 n1*p2*n3 ...
             n1*n2*p3 p1*n2*p3 p1*p2*p3 n1*p2*p3];
% Matrix of the gradient of shape functions 
BeXi = zeros(3,8);
BeXi(:, 1:4) = [-n2*n3  n2*n3  p2*n3 -p2*n3; ...
                -n1*n3 -p1*n3  p1*n3  n1*n3; ...
                -n1*n2 -p1*n2 -p1*p2 -n1*p2];
BeXi(:, 5:8) = [-n2*p3  n2*p3  p2*p3 -p2*p3; ...
                -n1*p3 -p1*p3  p1*p3  n1*p3; ...
                 n1*n2  p1*n2  p1*p2  n1*p2];
BeXi = BeXi * 1/8;
            