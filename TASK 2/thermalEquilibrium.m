function [fluxElement] = thermalEquilibrium(qheatGLO,CN,COOR)
%% medium heat at the element 
[nElement nPerElement ]= size(CN);
fluxElement = zeros(nElement,1);
for e = 1:nElement
    xDirec = -qheatGLO(1,e)-qheatGLO(3,e)+qheatGLO(5,e)+qheatGLO(7,e);
    yDirec = +qheatGLO(2,e)-qheatGLO(4,e)-qheatGLO(6,e)+qheatGLO(8,e);
    fluxElement(e) = abs(xDirec)+abs(yDirec);
end 
