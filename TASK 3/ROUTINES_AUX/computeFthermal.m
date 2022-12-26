function Fthermal = computeFthermal(COOR,CN,TypeElement, celasglo,alfa,tempNode);
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
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'K'; %Para octaedros esto no hace falta
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% Assembly of matrix K
% ----------------
Fthermal = zeros(nnode*ndim,1) ;

for e = 1:nelem
    celas = celasglo(:,:,e) ;  % Stiffness matrix of element "e"
    CNloc = CN(e,:) ;   % Coordinates of the nodes of element "e"
    Xe = COOR(CNloc,:)' ;     % Computation of elemental stiffness matrix
    beta = celas * alfa;
    tempNodeElement = tempNode(CNloc);
    FeThermal = ComputeFeThermalMatrix(tempNodeElement,weig,shapef,dershapef,Xe,beta);
   	%degress = Nod2DOF(CNloc,ndim);
    %Fthermal(degress, degress) =  Fthermal(degress, degress) + FeThermal;
    for anod=1:nnodeE 
            a = Nod2DOF(anod,ndim);  
            Anod = CN(e,anod) ; 
            A = Nod2DOF(Anod,ndim);     
            Fthermal(A) = Fthermal(A) + FeThermal(a); 
    end
    
    if mod(e,10)==0  % To display on the screen the number of element being assembled
        disp(['e=',num2str(e)])
    end
end