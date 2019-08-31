function T=tetrahedrization(Q)

% T=tetrahedrization(Q)
% 
% Input:
%    Q : basic quad partition
% Output:
%    T : tetrahedral partition with the same nodes and boundary 

T.coordinates=Q.coordinates;

% partition of each quad into 6 tetrahedra

tetra=[1 2 4 5;...
    6 4 5 8;...
    6 5 4 2;...
    2 3 4 6;...
    8 7 6 3;...
    3 8 4 6];
Ncubes=size(Q.elements,1);
T.elements=zeros(6*Ncubes,4);
for i=1:Ncubes
    elt=Q.elements(i,:);
    T.elements(6*(i-1)+(1:6),:)=elt(tetra);
end

% partition of each face into 2 triangles

tri=[1 2 4;2 3 4];

% Dirichlet faces

NDir=size(Q.dirichlet,1);
T.dirichlet=zeros(2*NDir,3);   
for i=1:NDir
    elt=Q.dirichlet(i,:);
    T.dirichlet(2*i-1:2*i,:)=elt(tri);
end

% Neumann faces

NNeu=size(Q.neumann,1);
T.neumann=zeros(2*NNeu,3);   
for i=1:NNeu
    elt=Q.neumann(i,:);
    T.neumann(2*i-1:2*i,:)=elt(tri);
end


return



