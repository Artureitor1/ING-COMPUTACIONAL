function reactionAdder(React,COOR)
%Reaction Forces
nnode = size(COOR,1);
ndim = size(COOR,2);
ReactOrdered = reshape(React,ndim,nnode);
Rx = sum(ReactOrdered(1,:));
Ry = sum(ReactOrdered(2,:));
Rz = sum(ReactOrdered(3,:));


%In order to locate the coordinates center at the centroid of the profile
%a tranfer of the points will be needed
auxiliaryCoord = [COOR';ones(1,nnode)];
transferMatrix = eye(ndim+1); % Transfer matrix to relocate the nodes
transferMatrix(:,end) = transferMatrix(:,end)+[0; -0.25/2; -0.25/2; 0];
coorTransfered = transferMatrix*auxiliaryCoord; %Coordinates relocated at the centroid of the profile

MomentsOrdered = zeros(ndim,nnode);

for node = 1:nnode
    d = coorTransfered([1,2,3],node)';
    F = ReactOrdered([1,2,3],node)';
    MomentsOrdered(:,node) = cross(d, F);
end

Mx = sum(MomentsOrdered(1,:));
My = sum(MomentsOrdered(2,:));
Mz = sum(MomentsOrdered(3,:));

disp(['Rx:', num2str(Rx)])
disp(['Ry:', num2str(Ry)])
disp(['Rz:', num2str(Rz)])

disp(['Mx:', num2str(Mx)])
disp(['My:', num2str(My)])
disp(['Mz:', num2str(Mz)])

end