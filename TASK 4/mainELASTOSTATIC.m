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
dampingFactor = 0.0;
 
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
t = linspace(0,5* 2*pi/Freq(5), 500);
DISP = zeros(length(d), length(t));

for i = 1:length(t)
    DISP(DOFr,i) = d(DOFr);
    DISP(DOFl,i) = dampedVibration(d, dampingFactor, MODES, Freq, M, DOFl, t(i));
end
%% 
%% Amplitude of modes graphic
qi0 = zeros(length(Freq), 1);
for i=1:length(Freq)
    qi0(i) = MODES(:,i)'*M(DOFl,DOFl)*d(DOFl);
end

figure();
bar(1:length(Freq), abs(qi0));

% PostPROCESS dynamic
GidPostProcessDynamic(COOR, CN, TypeElement, DISP, NAME_INPUT_DATA, posgp,...
    NameFileMesh, t)

%% Forced DAMPING
% dampingFactor = 0.99;
% t = linspace(0,40* 2*pi/Freq(5), 250);
% F_freq = linspace(1, 1000, 100);
% DISP_STEADY_MAX = zeros(length(F_freq), 1);
% DISP_MAX = zeros(length(F_freq), 1);
% for j = 1:length(F_freq)
%     DISP_forced = zeros(length(d), length(t));
%     DISP_steady = 0;
%     for i = 1:length(t)
%         DISP_forced(DOFr,i) = d(DOFr);
%         DISP_forced(DOFl,i) = dampedForcedVibration(d,dampingFactor,MODES, Freq, ...
%             DOFl, Ftrac(DOFl), F_freq(j), t(i));
% 
%         DISP_steady = max(DISP_steady, dampedForcedVibrationMaxDSteady(d,dampingFactor,MODES, Freq, ...
%             DOFl, Ftrac(DOFl), F_freq(j), t(i)));
%     end
%     DISP_MAX(j) = max(abs(DISP_forced),[],'all');
%     DISP_STEADY_MAX(j) = DISP_steady;
%     fprintf('iteration %d, frequency:%f\n', j, F_freq(j));
%     
% end
% 
% figure();
% hold on;
% plot(F_freq, 20*log10(DISP_STEADY_MAX./1e-12));
% plot(F_freq, 20*log10(DISP_MAX./1e-12))
% hold off;
% 
% %GidPostProcessDynamicDamped(DISP_forced,t,nodeID)


