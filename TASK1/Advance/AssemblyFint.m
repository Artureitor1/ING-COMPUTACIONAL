function [Residual,STRAIN,STRESS] = AssemblyFint(COOR,CN,d_k,stressFUN,AreaFUN) 

numberElement = size(CN,1) ; % Number of rows and columns SIZE
numberNode = numberElement+1 ;
nodesPerElement = size(CN,2) ;
Residual = zeros(numberNode,1);
STRAIN = zeros(numberElement,1);
STRESS = zeros(numberElement,1);
    for e=1:numberElement
        elementNodes = CN(e,:);    % Global numbering of nodes of element "e"
        elementCoors = COOR(elementNodes);
        he = elementCoors(2)-elementCoors(1);
        areaXcoord = AreaIntegral(he, elementCoors,AreaFUN); %Dominio gaus 
        bXcoord  = (1/he)*[-1 1];
        strainXcoord = stressFUN((1/he)*[-1 1]*d_k(CN(e,:)));
        stressXcoord = bXcoord*d_k(CN(e,:));
        FeXcoord = bXcoord'*areaXcoord*strainXcoord;
        Residual(CN(e,:),1) = Residual(CN(e,:),1)+FeXcoord;
        STRAIN(e,1) = strainXcoord';
        STRESS(e,1) = stressXcoord';
    end
end

