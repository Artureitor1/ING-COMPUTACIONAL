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

%% Response of one node with different damping ratios
% This code could be more optimized so it would only compute the
% displacement of one of the nodes instead of every node
dampingFactors = [ 0.05 0.1 0.3 0.5 0.99];
axis = 3; % x = 1, y = 2, z = 3
prevalent_mode = 5;
node = 38;
t = linspace(0,5* 2*pi/Freq(prevalent_mode), 150);
figure();
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times
hold on;
for j = 1:length(dampingFactors)
    DISP = zeros(length(d), length(t));
    for i = 1:length(t)
        DISP(DOFr,i) = d(DOFr);
        DISP(DOFl,i) = dampedVibration(d, dampingFactors(j), MODES, Freq, M, DOFl, t(i));
    end
    movement = DISP((38-1)*3+axis, :)';
    plot(t, movement, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\bar{\\xi}_i$ = %0.2f',dampingFactors(j)));
end
title('Displacement in free damped vibration - Bending Load');
xlabel('Time [s]') 
ylabel('Displacement [m]')
leg = legend();
set(leg,'interpreter','latex');
legend;
hold off
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


