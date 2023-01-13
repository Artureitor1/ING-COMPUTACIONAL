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
dampingFactor = 0.01;
 
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
[d, strainGLO, stressGLO, React, posgp, MODES, DOFl, M, Freq, Ftrac]= SolveElastFE(COOR,CN,TypeElement,TypeElementB, celasglo, ...
    densglo, DOFr,dR, Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA)  ; 

reactionAdder(React,COOR);

% POSTPROCESS
% --------------------------------------------
GidPostProcessModes(COOR,CN,TypeElement,MODES,posgp, NameFileMesh,DATA,DOFl)

GidPostProcess(COOR,CN,TypeElement,d,strainGLO, stressGLO,  React,NAME_INPUT_DATA,posgp,NameFileMesh,DATA);


% %% DAMPING VIBRATION
% t = linspace(0,40* 2*pi/Freq(5), 250);
% DISP = zeros(length(d), length(t));
% 
% for i = 1:length(t)
%     DISP(DOFr,i) = d(DOFr);
%     DISP(DOFl,i) = dampedVibration(d, dampingFactor, MODES, Freq, M, DOFl, t(i));
% end
% 
% % PostPROCESS dynamimc
% GidPostProcessDynamic(COOR, CN, TypeElement, DISP, NAME_INPUT_DATA, posgp,...
%     NameFileMesh, t)

%% Forced DAMPING
t = linspace(0,40* 2*pi/Freq(5), 250);
DISP_forced = zeros(length(d), length(t));
force_w = 5; % rad/s
for i = 1:length(t)
    DISP_forced(DOFr,i) = d(DOFr);
    DISP_forced(DOFl,i) = dampedForcedVibration(d,dampingFactor,MODES, Freq, ...
        DOFl, Ftrac(DOFl), force_w, t(i));
end

GidPostProcessDynamic(COOR, CN, TypeElement, DISP_forced, NAME_INPUT_DATA, posgp,...
    NameFileMesh, t)


