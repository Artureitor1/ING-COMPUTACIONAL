function [weig,posgp,shapef,dershapef] = Triangle3NInPoints(TypeIntegrand)
% This function returns, for each 2-node 1D linear element,
% and using a 1 points Gauss rule (ngaus=1): , if TypeIntegrand = K
% and ngaus=2, if TypeIntegrand = RHS:
%
% weig = Vector of Gauss weights (1xngaus)
% posgp: Position of Gauss points  (ndim x ngaus)
% shapef: Array of shape functions (ngaus x nnodeE)
% dershape: Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)

switch TypeIntegrand
    case {'K'}
        % One integration point
        weig  = [.5] ;
        posgp = [[1/3 1/3]];
    case {'RHS'}
        % Three integration points
        weig  = [1 1 1]./6 ;
        posgp = 0.5*[1 0; 1 1; 0 1] ;
end

ndim = 3; nnodeE = 3 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig)
    xi = posgp(g,:);
    [Ne, BeXi] = Triangle3N(xi) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end