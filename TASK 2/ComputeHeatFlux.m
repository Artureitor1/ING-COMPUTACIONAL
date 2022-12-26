function [qheatGLO, posgp]= ComputeHeatFlux(COOR,CN,TypeElement,ConductMglo,d) 

nnode = size(COOR,1);  % Number of nodes
ndim = size(COOR,2);   % Spatial Dimension of the problem  (2 or 3)
nelem = size(CN,1);   % Number of elements 
nnodeE = size(CN,2) ; %Number of nodes per element 
qheatGLO = zeros(nnodeE*ndim,nelem); 
TypeIntegrand = 'Quadrilateral';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;

for e = 1:nelem
    elementNodes = CN(e,:);
    elementCoord = COOR(elementNodes,:)';
    qE = ComputeHeatFluxMatrix(ConductMglo(:,:,e), weig, dershapef, elementCoord,d(CN(e,:)));
    qheatGLO(:,e) =  qE; 
end
end