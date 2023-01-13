function[d strainGLO stressGLO  React posgp MODES DOFl Mm FREQ]  = SolveElastFE(COOR,CN,TypeElement,TypeElementB, celasglo, densglo,  DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA) ; 

%%% This function returns the (nnode*ndim x 1) vector of nodal displacements (d),
%%% as well as the arrays containing  the stresses (stressGLO) and strains (strainGLO) at all gauss
%%% points 
% % INPUTS
% --------------
% 1. Finite element mesh 
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% TypeElementB: Type of boundary finite element (linear...)
% -----------
% 2. Material
% -----------
%  celasglo (nstrain x nstrain x nelem)  % Array of elasticity matrices
%  celasgloINV (6 x 6 x nelem)  % Array of compliance matrices (3D)
% -------------------------
% 3. Dirichlet (Essential) Boundary Condition s
% --------------------------------------------
%  DOFr --> Set of Global Degrees of Freedom with prescribed displacements 
%  dR   --> Vector of prescribed displacements  (size(DOFr) = size(dR))
% ---------------------------------------
% 4. Neumann (natural) Boundary Conditions
% -------------------------------------------
% DISTRIBUTED LOADS
% -------------------
%  CNb: Cell array in which the cell CNb{idim} contains the connectivity matrix for the boundary elements
%  of the traction boundaries in the idim direction
%  Tnod: Cell array in which the entry Tnod{idim} features the vector with the prescribed traction at
%   the nodes specified in CNb{idim}    (note that size(CNb{idim}) = size(Tnod{idim}))
%  each cell of 
%  POINT LOADS
% --------------------------------
%  Fpnt  (nnode*ndime x 1):  Vector containing point forces applied on the
%  nodes of the discretization
% ----------------------------------
% 5. Body force
% ---------------
%  fNOD: Vector containing the nodal values of the heat source function (nnode*ndime x1 )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=[]; strainGLO=[] ; stressGLO=[] ;posgp=[] ; 

neig = 20;
nnode = size(COOR,1); 
ndim = size(COOR,2); 
DOFl = 1:nnode*ndim;
DOFl(DOFr) = [] ;


% A) Global stiffness matrix 
% ------------------------------
disp('Computing stiffness matrix K ...')
K = ComputeK(COOR,CN,TypeElement, celasglo) ; 

disp('Computing mass matrix M ...')
Mm = ComputeM(COOR,CN,TypeElement, densglo) ; 

% B) External force vector due to body forces
% ------------------------------
disp('Computing   external force vector due to body forces (Fb)...')
Fb = ComputeFb(COOR,CN,TypeElement, fNOD);

% C)  External force vector due to   boundary tractions  
% ------------------------------
disp('Computing  external force vector due to   boundary tractions ..')
Ftrac = FtracCOMP(COOR,CNb,TypeElementB,Fpnt,Tnod);

% D)  Modal frequencies analysis   
% ------------------------------
disp('Computing  modal analysis...')
[MODES FREQ] = UndampedFREQ(Mm(DOFl,DOFl),K(DOFl,DOFl),20) 

% % E)  First PostProcesing  
% % ------------------------------
% GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,'Beam_4',DATA,DOFl)

% D) Solving for the vector of unknown displacements 
disp('Solving...')

[d strainGLO stressGLO  React posgp] = SolveELAS(K,Fb,Ftrac,dR,DOFr,COOR, ...
    CN,TypeElement,celasglo,densglo, typePROBLEM,celasgloINV,DATA) ; 


save('INFO_FE.mat','React','COOR','DOFr')