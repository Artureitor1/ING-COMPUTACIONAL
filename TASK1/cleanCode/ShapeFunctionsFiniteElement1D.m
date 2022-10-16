function [N, xNODES] = ShapeFunctionsFiniteElement1D(n,x0,xL)
syms x
if nargin ==0
    n = 5 ;
    x0 = 0;
    xL= 1;
end
nnodes = n+1;
xNODES = linspace(x0,xL,nnodes);
N = sym(zeros(1,nnodes)) ;
% figure(30)
% hold on
% xlabel('x')
% ylabel('N(x)')

for inode = 1:nnodes
    if inode ==1
        x1 = xNODES(inode) ; x2 = xNODES(inode+1) ;
        N(inode) = piecewise( (x1 <=x) & (x<=x2), (x2-x)/(x2-x1),0) ;
    elseif inode == nnodes
        xnp1 = xNODES(inode) ; xn = xNODES(inode-1) ;
        N(inode) = piecewise( (xn <=x) & (x<=xnp1), (x-xn)/(xnp1-xn),0) ;
    else
        xi = xNODES(inode) ; xim1 = xNODES(inode-1) ;  xip1 =  xNODES(inode+1) ; 
        N(inode) = piecewise( (xim1 <=x) & (x<=xi), (x-xim1)/(xi-xim1), (xi <=x) & (x<=xip1), (xip1-x)/(xip1-xi)  , 0);
    end
end


