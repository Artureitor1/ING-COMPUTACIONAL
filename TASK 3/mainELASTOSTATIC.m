clc
clear all
tic;
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
[COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,DATA] = ReadInputDataFile(NAME_INPUT_DATA)  ; 
% TEMPERATURE
tempNode = ones(length(COOR),1)*0;
alfa = [21;21;21;0;0;0]*10^(-6);
% SOLVER 
% --------------------------------------------
[d strainGLO stressGLO  React posgp]= SolveElastFE(COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA,tempNode,alfa)  ; 

reactionAdder(React,COOR);
%posgp = posgp';
% POSTPROCESS
% --------------------------------------------
GidPostProcess(COOR,CN,TypeElement,d,strainGLO, stressGLO,  React,NAME_INPUT_DATA,posgp,NameFileMesh,DATA);
toc;