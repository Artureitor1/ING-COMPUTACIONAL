function Mm = ComputeM(COOR,CN,TypeElement, densglo) ;
%%%%
% This subroutine   returns the global stiffness matrix K (ndim*nnode x ndim*nnode)
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE), % TypeElement: Type of finite element (quadrilateral,...),  celasglo (nstrain x nstrain x nelem)  % Array of elasticity matrices
% Dimensions of the problem
if nargin == 0
    load('tmp1.mat')
end
nnode = size(COOR,1); 
ndim = size(COOR,2); 
nelem = size(CN,1);
nnodeE = size(CN,2) ;  
vecDeg = (1:1:nnode*ndim); %vecDeg and matDeg are used to relate local to global DEG
matDeg = reshape(vecDeg,ndim,nnode);
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'K';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% Assembly of matrix K
% ----------------
Mm = sparse([],[],[],nnode*ndim,nnode*ndim,nnodeE*ndim*nelem);
for e = 1:nelem
    dens = densglo(e) ;  % Stiffness matrix of element "e"
    CNloc = CN(e,:) ;   % Coordinates of the nodes of element "e"
    Xe = COOR(CNloc,:)' ;     % Computation of elemental stiffness matrix
    Me = ComputeMeMatrix(dens,weig, shapef, dershapef, Xe) ;

    degress = reshape(matDeg(:,CNloc),nnodeE*ndim,1); %For didactical purposes it has been decided to use a local to global degree converter of our own (without using Nod2DOF).
    Mm(degress, degress) =  Mm(degress, degress) + Me;

    if mod(e,10)==0  % To display on the screen the number of element being assembled
        disp(['e=',num2str(e)])
    end
end


