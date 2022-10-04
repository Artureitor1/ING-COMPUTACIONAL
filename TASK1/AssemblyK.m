function K = AssemblyK(COOR, CN, rho)
    % COOR: Coordinate matrix % CN %Element connectivy matrix
    nelem = size (CN, 1);nnode = size (COOR, 1); ;nnodeE = size (CN, 2);
    K =sparse(nnode, nnode);
    for e=1:nelem
        % Element matrix
        NODOSe = CN(e,:);
        COOR_e = COOR (NODOSe);
        he= COOR_e(2)-COOR_e(1);
        % Elemental matrix
        kep = 
        Ke= 1/he*[1 -1; -1 1]-kep;
        % Assembly
        for a = 1:nnodeE
            for b = 1:nnodeE
                A = CN(e,a);
                B = CN(e,b);
                K(A,B)= K(A,B) + Ke(a,b); 
            end
        end
    end
end