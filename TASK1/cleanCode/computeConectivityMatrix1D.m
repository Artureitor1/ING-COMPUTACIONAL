function conectivityMatrix = computeConectivityMatrix1D(nodes)
    conectivityMatrix = zeros(length(nodes)-1,2);
    for n = 1:(length(nodes)-1)
        conectivityMatrix(n,1) = n;
        conectivityMatrix(n,2) =n+1;
    end 
end