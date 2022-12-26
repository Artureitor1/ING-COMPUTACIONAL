function [qE] = ComputeHeatFluxMatrix(ConductM,weig,dershapef,Xe,de)

ndim = size(Xe,1); 
ngaus = length(weig); 
nnodeE = size(Xe,2); 
qE = zeros(nnodeE,ndim); 

for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ; 
    % Jacobian Matrix 
    Je = Xe*BeXi' ; 
    % Matrix of derivatives with respect to physical coordinates 
    Be = inv(Je)'*BeXi ; 
    %
    qE(g,:) = -(ConductM*Be*de)'; 
end
qE = reshape(qE',[nnodeE*ndim,1]);

end