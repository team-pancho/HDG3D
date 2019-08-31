% This Script is based on cube(a,b,c) and used to create a Corner with
% Dirichlet and Neumann BCs:
% The following domain is (0,1)^3/(0.5,1)^3. The Neumann boundary are
% z=0, z=1 and z=0.5 (at the corner), y=0.
% March 25, 2013
clear
clc
%%
Q=cube(2,2,2);
Q.coordinates(27,:)=[];
Q.elements(8,:)=[];
Q.dirichlet(4,:)=[14 15 24 23];
Q.dirichlet(12,:)=[14 23  26 17];
Q.neumann(8,:)=[14 17 18 15];
%%
T1=tetrahedrization(Q);
Q=partition(Q);
T2=tetrahedrization(Q);
Q=partition(Q);
T3=tetrahedrization(Q);
Q=partition(Q);
T4=tetrahedrization(Q);
%%
T1=HDGgrid3d(T1,-1);
T2=HDGgrid3d(T2,-1);
T3=HDGgrid3d(T3,-1);
T4=HDGgrid3d(T4,-1);