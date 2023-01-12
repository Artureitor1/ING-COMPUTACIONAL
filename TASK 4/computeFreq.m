function [MODES FREQ]= computeFreq(COOR,CN,DOFr,)
neig = 20;
nnode = size(COOR,1); 
ndim = size(COOR,2); 
nelem = size(CN,1); 
nnodeE = size(CN,2); 
DOFl = 1:nnode*ndim;
DOFl(DOFr) = [] ;

[MODES FREQ] = UndampedFREQ(M(DOFl,DOFl),K(DOFl,DOFl),neig) 

end