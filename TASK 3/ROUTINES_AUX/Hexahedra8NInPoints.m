function [weig,posgp,shapef,dershapef] = Hexahedra8NInPoints()
% This function returns, for each 2-node 1D linear element, % Esta
% definicion no esta bien 
% and using a 1 points Gauss rule (ngaus=1): , if TypeIntegrand = K
% and ngaus=2, if TypeIntegrand = RHS:
%
% weig = Vector of Gauss weights (1xngaus)
% posgp: Position of Gauss points  (ndim x ngaus)
% shapef: Array of shape functions (ngaus x nnodeE)
% dershape: Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)

% weig  = [1, 1, 1, 1, 1, 1, 1, 1] ;
% posgp = 1/sqrt(3)*[-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; ...
%                    -1 -1  1; 1 -1  1; 1 1  1; -1 1  1] ;
% 
% 
% ndim = 3; 
% nnodeE = 8;
% ngaus = length(weig) ;
% shapef = zeros(ngaus,nnodeE) ;
% dershapef = zeros(ndim,nnodeE,ngaus) ;
% for g=1:length(weig)
%     xi = posgp(g,:) ;
%     [Ne, BeXi] = Hexahedra8N(xi) ;
%     shapef(g,:) = Ne ;
%     dershapef(:,:,g) = BeXi ;
% end


ngaus = 8;
ndim = 3;
nnodeE = 8;

signs = [-1 -1 -1;
          1 -1 -1;
          1  1 -1;
         -1  1 -1;
         -1 -1  1;
          1 -1  1;
          1  1  1;
         -1  1  1].';
     
posgp = 1/sqrt(3)*signs;
weig = ones(1,8);

shapef = 1/(2^ndim)*ones(ngaus,nnodeE);
        
for i = 1:ngaus
    for j = 1:nnodeE
        for k = 1:ndim
            shapef(i,j) = shapef(i,j)*(1+signs(k,j)*posgp(k,i));
        end
    end
end

dershapef = zeros(ndim,nnodeE,ngaus);

for i = 1:nnodeE
   dershapef(:,:,i) = 1/(2^ndim)*signs;
   
   for row = 1:ndim
       for col = 1:nnodeE
           range = 1:ndim;
           for j = range(range ~= row)
               dershapef(row, col,i) = dershapef(row, col,i)*(1+signs(j,col)*posgp(j,i));
           end
       end
   end
end
                               
end