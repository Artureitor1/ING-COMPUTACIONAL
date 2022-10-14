function Ff = AssemblyF(COOR,CN, f, N)

nelem = size(CN,1) ; % Number of rows and columns SIZE
nnode = nelem+1 ;
nnodeE = size(CN,2) ;

Ff = zeros ( nnode , 1 );
    for e=1: nelem
        NODOSe = CN(e,:);    % Global numbering of nodes of element "e"
        COOR_e = COOR(NODOSe);
        he = COOR_e(2)-COOR_e(1);
        Fe = (he /2 )* COMPUTE_Fe_FORCE (f,N, COOR_e);
        for a=1: nnodeE
            A = CN( e , a );
            Ff (A) = Ff (A) + Fe ( a) ;
        end
    end
end
