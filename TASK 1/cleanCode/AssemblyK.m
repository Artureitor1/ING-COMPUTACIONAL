function stiffnesMatrix = AssemblyK(coor,conectivityMatrix, rho)
% COOR: Coordinate matrix
% CN % Element connectivy matrix  ;
numberElements = size(conectivityMatrix,1) ; % Number of rows and columns SIZE
numberNodes = conectivityMatrix(end,2);
nodePerElement = size(conectivityMatrix,2) ;
stiffnesMatrix =sparse(numberNodes, numberNodes);
for element=1:numberElements  % Loop over number of elements 
    % Element matrix
    elementNodes = conectivityMatrix(element,:);    % Global numbering of nodes of element "e"
    elementCoor = coor(elementNodes) ;
    he = elementCoor(2)-elementCoor(1) ; % Size finite element
    % Elemental matrix
    Ke = 1/he*[1 -1; -1 1] - rho * he/2 * [2/3 1/3; 1/3 2/3];
    % Assembly
    for a = 1:nodePerElement
        for b = 1:nodePerElement
            A = conectivityMatrix(element,a);
            B = conectivityMatrix(element,b);
            stiffnesMatrix(A,B) = stiffnesMatrix(A,B) + Ke(a,b);
        end
    end
    
end