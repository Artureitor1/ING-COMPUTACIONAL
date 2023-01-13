clc
clear all
% Finite Element Program for Elastostatic problems  
% ECA.
% Technical University of Catalonia
% JoaquIn A. Hdez, October 23-th, 2015
% ---------------------------------------------------
if exist('ElemBnd')==0
    addpath('ROUTINES_AUX') ,
end

 
%%% INPUT  %%% 
% Input data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NAME_INPUT_DATA = 'BEAM3D' ;  % Name of the mesh file 
%NAME_INPUT_DATA_Termal = 'Thermal3D' ;  % Name of the mesh file 

%------------------------------------------------------

% PREPROCESS
[COOR,CN,TypeElement,TypeElementB, celasglo, densglo ,DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,DATA] = ReadInputDataFile(NAME_INPUT_DATA)  ; 


% SOLVER 
% --------------------------------------------
[d, strainGLO, stressGLO, React, posgp]= SolveElastFE(COOR,CN,TypeElement,TypeElementB, celasglo, ...
    densglo, DOFr,dR, Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA,tempNode,alfa)  ; 

reactionAdder(React,COOR);

% POSTPROCESS
% --------------------------------------------
GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp,'Beam_4',DATA,DOFl)