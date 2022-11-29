function [Residual,STRAIN,STRESS] = AssemblyFint(COOR,CN,d_k,stressFUN,AreaFUN)
numberElement = size(CN,1) ; % Number of rows and columns SIZE
numberNode = numberElement+1 ;
nodesPerElement = size(conectivityMatrix,2) ;
    for i = 1:numberElement
        elementNode = conectivityMatrix(e,:);    % Global numbering of nodes of element "e"
        elementCoor = coor(elementNode);
        he = elementCoor(2)-elementCoor(1);
        Fe = (he /2)* COMPUTE_Fe_FORCE(f,nodes, elementCoor);
        for a=1: nodesPerElement
            A = conectivityMatrix(e ,a);
            Ff (A) = Ff (A) + Fe (a) ;
        end
    end 
end 
