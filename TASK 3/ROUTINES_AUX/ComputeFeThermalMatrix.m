function FeThermal = ComputeFeThermalMatrix(tempNode,weig,shapef,dershapef,Xe,beta) ;
% Given 
% fe: Nodal values of the body force    (nnodeE*ndim x1)
% weig :   Vector of Gauss weights (1xngaus)
% shapef:   Array with the   shape functions at each Gauss point (ngaus x nnodeE )
% dershapef:   Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)
% Xe: Global coordinates of the nodes of the element,  
% % this function returns the element body force vector  Fbe
ndim = size(Xe,1) ; ngaus = length(weig) ; nnodeE = size(Xe,2)  ; 
FeThermal = zeros(ndim*nnodeE,1) ; 
for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ; 
    % Matrix of shape functions at point "g"
    NeSCL = shapef(g,:) ; 
    % Jacobian Matrix 
    Je = Xe*BeXi' ; 
    % JAcobian 
    detJe = det(Je) ;    
    %
    % Matrix of derivatives with respect to physical coordinates 
    BeTILDE = inv(Je)'*BeXi ; 
    % Matrix of derivatives with respect to physical coordinates 
    Be = QtransfB(BeTILDE,ndim) ;
    %
    FeThermal = FeThermal + weig(g)*detJe*Be'*(beta*(NeSCL*tempNode)); 
end