function [COOR,CN,TypeElement,CONNECTb,TypeElementB, ConductMglo,...
    dR,qFLUXglo,CNb,fNOD,NameFileMesh] = ...
    PreProcessInputData(NameFileMeshDATA,PROPMAT,DIRICHLET,NEUMANN,...
    fSOURCE)



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
% 2. MATERIAL PROPERTIES: output ConductMglo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConductMglo = zeros(ndim,ndim,nelem) ; 
% Conductivity matrix (isotropic)
for imat = 1:length(PROPMAT)
    kappa = PROPMAT(imat).kappa ;
    ConductM = kappa*eye(ndim) ; % eye = IDENTITY MATRIx
    ELEMS = find(MaterialType == imat) ;
    for eLOC=1:length(ELEMS)
        e = ELEMS(eLOC) ;
        ConductMglo(:,:,e) = ConductM ;
    end
end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of nodes at which temperature is prescribed
%  Specify the number of line(s) in which the temperature is imposed
rnod = cell(length(DIRICHLET),1)  ; dR =  cell(length(DIRICHLET),1) ; 
for  icond = 1:length(DIRICHLET)
    rnod{icond} =  ListOfNodesFACES(NameFileMesh,DIRICHLET(icond).NUMBER_SURFACE, ndim) ;
    dR{icond} = DIRICHLET(icond).PRESCRIBED_TEMPER*ones(size( rnod{icond} )) ; 
end
% Removed repeated condions 
% ---------------------------
rnod = cell2mat(rnod) ; 
dR = cell2mat(dR) ; 
[rnod, AAA] = unique(rnod) ;
dR = dR(AAA) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: qFLUXglo, CNb 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS
% ------------------------
CNb =cell(ndim,1) ; qFLUXglo=cell(ndim,1) ;
for  icond = 1:length(NEUMANN)
    rnodLOC=  ListOfNodesFACES(NameFileMesh,NEUMANN(icond).NUMBER_SURFACE,ndim) ;
    PRESCRIBED_qBAR = NEUMANN(icond).PRESCRIBED_qBAR;
    for idim = 1:ndim
        if PRESCRIBED_qBAR(idim) ~= 0
            t = PRESCRIBED_qBAR(idim) ;
            CNbLOC = ElemBnd(CONNECTb,rnodLOC) ;
            TnodLOC = t*ones(size(CNbLOC)) ;
            CNb{idim} = [CNb{idim} ; CNbLOC] ;
            qFLUXglo{idim} = [qFLUXglo{idim}  ; TnodLOC] ;
        end
    end
end
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Heat source : OUTPUT : fNOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Global vector of heat sources (constant)
fNOD = fSOURCE*ones(nnode,1) ; 

