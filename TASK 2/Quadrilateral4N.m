function  [Ne BeXi] = Quadrilateral4N(xi); 
% Shape functions and derivatives for 2-node linear element 
 
% Matrix of shape functions
Ne =0.25*[(1-xi(1))*(1-xi(2)) (1+xi(1))*(1-xi(2)) (1+xi(1))*(1+xi(2)) (1-xi(1))*(1+xi(2))];
% Matrix of the gradient of shape functions 
BeXi = 0.25*[-(1-xi(2)) (1-xi(2)) (1+xi(2)) -(1+xi(2));-(1-xi(1)) -(1+xi(1)) (1+xi(1)) (1-xi(1))]; 
