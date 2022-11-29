clear all
close all
elementsNumberArray = [4, 6, 8, 10];
errorArray= zeros(size(elementsNumberArray,2),2); %[displacement error, derivative error]
exactDisplacement = @(x) (61740346358469708232516810123*cos((2778046668940015^(1/2)*x)/16777216))/6174034635846970641698934560180 + (8773830921664641*x^2)/88897493406080480 + (2778046668940015^(1/2)*sin((2778046668940015^(1/2)*x)/16777216)*(66293192113335022366330415812131684352*2778046668940015^(1/2)*sin(2778046668940015^(1/2)/16777216) - 10977124066582370799746477992578602748793214053))/(18416558152385082682136054268153845979438084432645324800*cos(2778046668940015^(1/2)/16777216)) - 154350865896174268311882694656/7717543294808713302123668200225;
exactDerivative = @(x) (8773830921664641*x)/44448746703040240 - (61740346358469708232516810123*2778046668940015^(1/2)*sin((2778046668940015^(1/2)*x)/16777216))/103583112677085969401441632086004858880 + (cos((2778046668940015^(1/2)*x)/16777216)*(66293192113335022366330415812131684352*2778046668940015^(1/2)*sin(2778046668940015^(1/2)/16777216) - 10977124066582370799746477992578602748793214053))/(111221520341491811789912126265563782046673797120*cos(2778046668940015^(1/2)/16777216));

for n = 1:size(elementsNumberArray,2)
    elementNumber= elementsNumberArray(n);
    [nodeDisplacement,conectivityMatrix,coords,L]  = computeFEMmethod1D(elementNumber);
    errorArray(n,:) = computeError(nodeDisplacement,exactDisplacement,exactDerivative,coords,elementNumber,conectivityMatrix,L); 
end

% Plot

figure
hold on
title('Error of displacement and derivative', 'Interpreter','latex', 'FontSize',16)
plot(log(1./elementsNumberArray),log(errorArray(:,1)),LineWidth=2)
plot(log(1./elementsNumberArray),log(errorArray(:,2)),LineWidth=2)

grid on
grid minor
xlabel('Number of Elements', 'Interpreter','latex', 'FontSize',15)
ylabel('Error', 'Interpreter','latex', 'FontSize',15)
legend('Displacement', 'Derivative', 'location', 'northeast', 'Interpreter', 'latex', 'FontSize',12)
hold off


%Functions
function [nodeDisplacement,conectivityMatrix,coords,L] = computeFEMmethod1D (elementsNumber)
    %% Data
    %Problem Data
    syms x
    L = 1;
    g = 0.01;
    rho = (pi/L)^2;
    s = g*rho^2;
    f =  -s*x^2;
    F_AE = (g*pi^2)/L;
    
    %Numeric Data
    [nodes, coords] = ShapeFunctionsFiniteElement1D(elementsNumber, 0, L);
    nodes = 1:length (coords);
    restringedNodes = [1 -g];
    restringedForce= [nodes(end) F_AE];
    conectivityMatrix = computeConectivityMatrix1D(nodes);
    %Solver
    nodeDisplacement = computeDisplacement(nodes,restringedNodes,restringedForce,coords,conectivityMatrix,f,F_AE,rho);
end
function error = computeError(nodeDisplacement,exactDisplacement,exactDerivative,coords,elementNumber,conectivityMatrix,L)
    xi = [-1/sqrt(3) 1/sqrt(3)];
    ws = [1 1];
    error = [0 0];
    for element = 1:elementNumber
            errorPerElement = [0 0];
            elementDisplacement = [nodeDisplacement(conectivityMatrix(element,1));nodeDisplacement(conectivityMatrix(element,2))];
            elementCoords = [coords(conectivityMatrix(element,1)), coords(conectivityMatrix(element,2))];
            elementLength = elementCoords(2)-elementCoords(1);
            for g=1:2
                N = 1/2 * [1-xi(g) 1+xi(g)];
                B = 1/elementLength *[-1 1];
                integratingErrorDisplacement = (exactDisplacement(N * elementCoords')-(1/2 * [1-xi(g) 1+xi(g)])*elementDisplacement)^2;
                integratingErrorDerivative = (exactDerivative(N * elementCoords')-(B*elementDisplacement))^2;
                errorPerElement(1,1) = errorPerElement(1,1) + ws(g) * (integratingErrorDisplacement);
                errorPerElement(1,2) = errorPerElement(1,2) + ws(g) * (integratingErrorDerivative);
            end
            errorPerElement = (elementLength/2)*errorPerElement;
            error = error + errorPerElement;
    end
    error = sqrt(error);
end 
