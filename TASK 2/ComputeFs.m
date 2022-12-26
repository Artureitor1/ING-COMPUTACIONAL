function Fs = ComputeFs(COOR,CN,TypeElement, fNOD) ; 
% This subroutine   returns the  heat source contribution (Fs)    to the
% global flux vector. Inputs
% --------------
% 1. Finite element mesh 
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% -----------
% 2. Vector containing the values of the heat source function at the nodes
% of the mesh
% -----------
%  fNOD (nnode x 1)  %  
%%%%
 
% Dimensions of the problem 
nnode = size(COOR,1);
ndim = size(COOR,2);
nelem = size(CN,1);
nnodeE = size(CN,2) ;    

[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,'RHS') ; 

  
Fs = zeros(nnode,1) ;  
%......

for e = 1:nelem
    elementNodes = CN(e,:);
    elementCoord = COOR(elementNodes,:)';
    Fse = ComputeFseVector(fNOD(e), weig, shapef, dershapef, elementCoord);
    for j = 1:size(elementNodes)
        A = CN(e,j);
        Fs(A) =  Fs(A) + Fse(j);
    end
end