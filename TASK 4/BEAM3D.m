% Inputs example assigment 3

% --------------------------------------
% % 1. NAME OF THE MESH AND DATA FILES. COORDINATES, CONNECTIVITIES,  LISTS
% % OF NODES FOR IMPOSING BOUNDARY CONDITIONS
% %--------------------------------------------------------------------------
NameFileMeshDATA = 'plane' ;   %  This is the name of GIDs project
% wherein the mesh has been constructed. To generate the files needed by matlab,
% remember to export both the mesh and the conditions data. (Files >
% Export > Gid Mesh) and (Files > Export > Calculation File)
% ---------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
% 2. Type of structural problem (plane stress ('pstress'), plane strain ('pstrain'), 3D)
% --------------------------------------------------------------
typePROBLEM ='3D' ;
% -----------------------------------------------------------------------------------


% -----------------------------------------------------------------------------------
% 3. Material data
% -----------------------------------------------------------------------------------
imat =1 ; % Index material
% Elasticity matrix
E = 70000  ; %  MPa, Young's modulus
nu = 0.3; % Poisson's coefficient
% Compliance matrix for an isotropic materials (with all entries, 3D)
% See slides, page 23.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = E/2/(1+nu) ;  % Shear modulus
celasINV3D = [1/E  -nu/E -nu/E  0 0 0
    -nu/E  1/E  -nu/E  0 0 0
    -nu/E  -nu/E  1/E  0 0 0
    0      0    0  1/G   0 0
    0      0      0  0 1/G 0
    0      0      0  0  0  1/G] ;
ElasticityMatrix = inv(celasINV3D) ;
PROPMAT(imat).ElasticityMatrix =  ElasticityMatrix  ; %
% imat = 2 .... Repeat the above sequence of operations for defining
% another material

% -------------------------
%% ---5.3)  Body forces
% ---------------------
fBODY = 0 ;  % Constant value per unit volum MN/m^3.


% DENSITY (THIS IS FOR ASSIGNMENT 4)
dens0 = 2.7/1000;   %

%
DATA.PRINT_AVERAGE_STRESSES_ON_ELEMENTS = 0  ; % Print volumetric average of stresses at each element




% PROCESSING INPUT DATA

[COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,celasglo,DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,densglo,celasgloINV] = ...
    PreProcessInputData(NameFileMeshDATA,PROPMAT,...
    fBODY,dens0,typePROBLEM);



