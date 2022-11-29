function K = AssemblyKnon(COOR,CN,d_k,AreaFUN,DerStressFUN)
numberElement = size(CN,1) ; % Number of rows and columns SIZE
numberNode = numberElement+1 ;
nodesPerElement = size(CN,2) ;
K = zeros(numberNode,numberNode);
    for e=1:numberElement
        elementNodes = CN(e,:);    % Global numbering of nodes of element "e"
        elementCoors = COOR(elementNodes);
        he = elementCoors(2)-elementCoors(1);
        areaXcoord = AreaIntegral(he, elementCoors,AreaFUN); %Dominio gaus 
        bXcoord  = (1/he)*[-1 1];
        strainDerXcoord = DerStressFUN(bXcoord*d_k(CN(e,:)))*bXcoord;
        KeXcoord = bXcoord'*areaXcoord*strainDerXcoord;
        K(CN(e,:),CN(e,:)) = K(CN(e,:),CN(e,:))+KeXcoord;
    end
end 