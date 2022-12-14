function [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints()
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

weig  = [1 1 1 1] ;
posgp = 1/sqrt(3)*[-1 -1; 1 -1; 1 1; -1 1] ;


ndim = 2; 
nnodeE = 4;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig);
    xi = posgp(g,:) ;
    [Ne, BeXi] = Quadrilateral4N(xi) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end