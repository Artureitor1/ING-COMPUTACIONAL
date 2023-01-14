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
axis = 2; % x = 1, y = 2, z = 3
prevalent_mode = 1;
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
    movement = DISP((node-1)*3+axis, :)';
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
dampingFactor = 0.03;
t = linspace(0,30* 2*pi/Freq(5), 500);
F_freqs = [Freq(5)*0.8, Freq(5)*.9];
figure();
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times
hold on;
for j = 1:length(F_freqs)
    DISP_forced = zeros(length(d), length(t));
    for i = 1:length(t)
        DISP_forced(DOFr,i) = d(DOFr);
        DISP_forced(DOFl,i) = dampedForcedVibration(d,dampingFactor,MODES, Freq, ...
            DOFl, Ftrac(DOFl), F_freqs(j), t(i));
    end
    plot(t, DISP_forced((node-1)*3 + axis,:), ...
        'DisplayName',sprintf('$\\bar{\\omega}_i$ = %0.2f rad/s',F_freqs(j)));
end

title('Forced vibration: $\bar{\xi}_i = 0.03$ - Torque Load', 'Interpreter', 'latex');
xlabel('Time [s]') 
ylabel('Displacement in Z [m]')
leg = legend();
set(leg,'interpreter','latex');
legend;
hold off;

GidPostProcessDynamic(COOR, CN, TypeElement, DISP_forced, NAME_INPUT_DATA, posgp,...
    'Beam_4_disp.msh', t)

%% Frequency analysis


dampingFactors = [0.01, 0.1, 0.5, 0.99];
F_freq = linspace(1, 20000, 500);
DISP_STEADY_MAX = zeros(length(F_freq), 1);
figure();
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times
hold on;
for i = 1:length(dampingFactors)
    for j = 1:length(F_freq)
        d1 = zeros(length(d),1);
        d1(DOFl) = dampedForcedVibrationMaxDSteady(d,dampingFactors(i),MODES, Freq, ...
                DOFl, Ftrac(DOFl), F_freq(j), t(i));
        DISP_STEADY_MAX(j) = abs(d1((node-1)*3 + axis));
    end
    plot(F_freq, 20*log10(DISP_STEADY_MAX./1e-12), ...
        'DisplayName',sprintf('$\\bar{\\xi}_i$ = %0.2f',dampingFactors(i)));
end
title('Frequency response of the beam - Torque Load');
xlabel('Frequency [rad/s]') 
ylabel('Displacement [dB]')
leg = legend();
set(leg,'interpreter','latex');
legend;
hold off;


