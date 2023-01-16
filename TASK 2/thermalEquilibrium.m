function [fluxElement] = thermalEquilibrium(qheatGLO,CN,COOR)
%% medium heat at the element 
[nElement nPerElement ]= size(CN);
fluxElement = zeros(nElement,1);
for e = 1:nElement
    xDirec = -qheatGLO(1,e)-qheatGLO(4,e)-qheatGLO(7,e)-qheatGLO(10,e)+qheatGLO(13,e)+qheatGLO(16,e)+qheatGLO(19,e)+qheatGLO(22,e);
    yDirec = -qheatGLO(2,e)+qheatGLO(5,e)+qheatGLO(8,e)-qheatGLO(11,e)-qheatGLO(14,e)+qheatGLO(17,e)+qheatGLO(20,e)-qheatGLO(23,e);
    zDirec = +qheatGLO(3,e)+qheatGLO(6,e)-qheatGLO(9,e)-qheatGLO(12,e)+qheatGLO(15,e)+qheatGLO(18,e)-qheatGLO(21,e)-qheatGLO(24,e);
    fluxElement(e) = abs(xDirec)+abs(yDirec)+abs(zDirec);
end 
