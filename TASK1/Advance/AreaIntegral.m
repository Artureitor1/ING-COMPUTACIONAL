function [Area] = AreaIntegral(he, elementCoor,AreaFUN)
%Quadrature Gaus integral method is used. Order 2 and 2 points  
xi = [-1/sqrt(3) 1/sqrt(3)];
ws = [1 1];

coordQuadrature = @(xi) 0.5*((1-xi)*elementCoor(1)+(1+xi)*elementCoor(2));
integralQuadrature  = 0;
for g=1:2
    AreaQuadrature = AreaFUN(coordQuadrature(xi(g)));
    integralQuadrature = AreaQuadrature*ws(g)+integralQuadrature;
end
Area = integralQuadrature*(he/2);
end