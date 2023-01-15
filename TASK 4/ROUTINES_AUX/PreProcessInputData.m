function [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,celasglo,...
    DOFr,dR,Tnod,CNb,fNOD,Fpnt,NameFileMesh,densglo,celasgloINV] = ...
    PreProcessInputData(NameFileMeshDATA,PROPMAT,...
    fBODY,dens0,typePROBLEM)



% ----------------------------
%%%%%%%%%%%%%%%%%
 % ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB
NameFileMesh = [NameFileMeshDATA,'.msh']; % Name of the file containing the mesh information (Generated with GID)
[COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
    ReadMeshFile(NameFileMesh)  ;

nnode = size(COOR,1) ;% Number of nodes 
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%typePROBLEM = 'pstress';  %'pstress'/'pstrain'/'3D';  Plane stress/ plane strain problem 
if ndim==2 
    nstrain = 3; 
else
    nstrain = 6 ; 
    typePROBLEM ='3D' ;
end
celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
celasgloINV = zeros(6,6,nelem) ;

for imat = 1:length(PROPMAT)    
    celas3D =PROPMAT(imat).ElasticityMatrix ; %
    INVcelas3D = inv(celas3D) ; 
    ELEMS = find(MaterialType == imat) ;
    
    switch typePROBLEM
        case 'pstrain'
            rowcol = [1 2 6] ;
            celas = celas3D(rowcol,rowcol) ;
        case 'pstress'
            rowcol = [1 2 6] ;
            celasINV3D = inv(celas3D) ;
            celasINV = celasINV3D(rowcol,rowcol) ;
            celas = inv(celasINV) ;
        case '3D'
            celas = celas3D ;
    end
    for eLOC=1:length(ELEMS)
        e = ELEMS(eLOC) ;
        celasglo(:,:,e) = celas ;
        celasgloINV(:,:,e) = INVcelas3D ; 
    end
end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of nodes at which displacement is prescribed (in any of the x-y and z directions)
rnod = cell(ndim,1) ; uPRES=cell(ndim,1)  ; 

 
% Degrees of freedom and prescribed displacements 
DOFr = [] ; dR = [] ; 
 
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS
% ------------------------
CNb =cell(ndim,1) ; Tnod=cell(ndim,1) ;
      

% POINT LOADS 
% -----------
Fpnt =zeros(ndim*nnode,1) ; % There are no point loads

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Body forces
%% 
fNOD = fBODY*ones(nnode*ndim,1) ; 

% 6. Density 
densglo = dens0*ones(nelem,1);  % in kKg/m^3 


% 
