clc
clear all
% "Template" Finite Element Program for Heat Conduction Problems
% ECA.
% Technical University of Catalonia
% JoaquIn A. Hdez, October 8-th, 2015
% ---------------------------------------------------
if exist('ElemBnd')==0
    addpath('ROUTINES_AUX') ;
end 
%%% INPUT  %%% 
% Input data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NAME_INPUT_DATA = 'DATA_ASSIGNMENT2';
%------------------------------------------------------

% PREPROCESS  
[COOR,CN,TypeElement,TypeElementB, ConductMglo,  rnod,dR,...  
    qFLUXglo,CNb,fNOD,NameFileMesh] = ReadInputDataFile(NAME_INPUT_DATA)  ; 

% SOLVER 
% --------------------------------------------
[d,qheatGLO,posgp] = SolveHeatFE(COOR,CN,TypeElement,TypeElementB, ConductMglo,  rnod,dR,...  
    qFLUXglo,CNb,fNOD)  ; 
[fluxElement] = thermalEquilibrium(qheatGLO,CN,COOR)

% POSTPROCESS
% --------------------------------------------
disp('POSTPROCESS....')
GidPostProcess(COOR,CN,TypeElement,d,qheatGLO,NAME_INPUT_DATA,posgp,NameFileMesh);
[elements, ~] = size(CN);
figure
hold on
title('Thermal equilibrium error per element', 'Interpreter','latex', 'FontSize',16)
bar(1:elements, fluxElement)
grid on
grid minor
xlabel('Element', 'Interpreter','latex', 'FontSize',15)
ylabel('HeatFlux (W/m2)', 'Interpreter','latex', 'FontSize',15)
hold off