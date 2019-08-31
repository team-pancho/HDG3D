function Q=cube(Nx,Ny,Nz)

% function Q=cube(Nx,Ny,Nz)
% input
%    Nx  : number of partitions in the x direction
%    Ny  : number of partitions in the y direction
%    Nz  : number of partitions in the z direction
% output
%    Q   : quad partition data structure
%          the top and bottom faces are Dirichlet

% see the program for possible modifications   

Nx=Nx+1;
Ny=Ny+1;
Nz=Nz+1;

% creating the cubic partition

list=reshape(1:Nx*Ny*Nz,Nx,Ny,Nz);
c=list(1:Nx-1,1:Ny-1,1:Nz-1);
q=[c(:) c(:)+1 c(:)+Nx+1 c(:)+Nx];
Q.elements=[q q+Nx*Ny];

% coordinates of the nodes
% you can change the cube to another paralellepiped here

x=linspace(0,1,Nx);        
y=linspace(0,1,Ny);
z=linspace(0,1,Nz);
[y,x,z]=meshgrid(y,x,z);
Q.coordinates=[x(:) y(:) z(:)];

reverse=[1 4 3 2];

% faces 1 (bottom z=0) & 2 (top z=1)
c=list(1:Nx-1,1:Ny-1,1);
face1=[c(:) c(:)+1 c(:)+Nx+1 c(:)+Nx];
face2=(Nz-1)*Nx*Ny+face1;
face2=face2(:,reverse);

% faces 3 (front y=0) & 4 (back y=1)

c=list(1:Nx-1,1,1:Nz-1);
face3=[c(:) c(:)+Nx*Ny c(:)+Nx*Ny+1 c(:)+1];
face4=(Ny-1)*Nx+face3;
face4=face4(:,reverse);

% faces 5 (left x=0) & 6 (right x=1)

c=list(1,1:Ny-1,1:Nz-1);
face5=[c(:) c(:)+Nx c(:)+Nx*Ny+Nx c(:)+Nx*Ny];
face6=(Nx-1)+face5;
face6=face6(:,reverse);

% The Dirichlet/Neumann division can be changed here

Q.dirichlet=[face4;face5;face6];%[face1;face2;face3;face4;face5;face6];
Q.neumann  =[face1;face2;face3];

