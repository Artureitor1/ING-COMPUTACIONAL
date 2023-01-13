function [MODES FREQ]= computeFreq(COOR,CN,DOFr,M,K)


[MODES FREQ] = UndampedFREQ(M(DOFl,DOFl),K(DOFl,DOFl),neig) 

end