% Inputs example assigment 2 
% ----------------------------
%--------------------------------------------------------------------------
% 1. NAME OF THE MESH AND DATA FILES. COORDINATES, CONNECTIVITIES,  LISTS
% OF NODES FOR IMPOSING BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
NameFileMeshDATA = 'Beam_4' ;   % For reading data from GID.
% -----------------------------------------------------------
% 2. Material data. 
% --------------------------------------------------------------
imat =1 ; % Index material
PROPMAT(imat).kappa =  10  ; %  Conductivity  of material "imat" (ISOTROPIC)

% -----------------------------------------------------------
% 3. Dirichlet boundary conditions (prescribed temperature)
% -----------------------------------------------------------
icond = 1; % Number of condition
DIRICHLET(icond).NUMBER_SURFACE = 1 ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
DIRICHLET(icond).PRESCRIBED_TEMPER = 100 ;

% -------------------------------------------------
% 4. Neumann Boundary conditions (prescribed flux)
% ------------------------------------------------
icond= 1 ;
NEUMANN(icond).NUMBER_SURFACE = 2 ;  % Surface on which the load is applied
NEUMANN(icond).PRESCRIBED_qBAR= [0.05 0 0] ;

% -------------------------------------------
% 4. Heat source (constant all over the body)
% --------------------------------------------
fSOURCE = 20; 


% END INPUTS 
% ----------------------------------------------------------------------------------------------

[COOR,CN,TypeElement,CONNECTb,TypeElementB, ConductMglo,...
    dR,rnod,qFLUXglo,CNb,fNOD,NameFileMesh] = ...
    PreProcessInputData(NameFileMeshDATA,PROPMAT,DIRICHLET,NEUMANN,...
    fSOURCE)


