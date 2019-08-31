function [uhd,qhn]=BC3d(uD,gx,gy,gz,T,k,formula)

%[uhd,qhn]=BC3d(uD,gx,gy,gz,T,k,formula)
%
% Input:
%          uD: Dirichlet data; vectorized function of three variables
%    gx,gy,gz: Neumann data (corresponds to kappa*grad(u))
%              vectorized functions of three variables
%           T: Full tetrahedrization
%           k:  polynomial degree
%     formula: quadrature formula in 2d (N x 4 matrix)
%
% Output:
%         uhd: d2 x Ndir,   Ndir, number of Dirichlet faces
%         qhn: d2 x Nneu,   Nneu, number of Neumann faces
%
% Last modified: March 21, 2013

x=T.coordinates(:,1);
y=T.coordinates(:,2);
z=T.coordinates(:,3);

%Dirichlet
xx=formula(:,[1 2 3])*x(T.dirichlet');  % Nnodes x Ndir
yy=formula(:,[1 2 3])*y(T.dirichlet');  % Nnodes x Ndir
zz=formula(:,[1 2 3])*z(T.dirichlet');  % Nnodes x Ndir
D=dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k);  % Nnodes x d2
wD=bsxfun(@times,formula(:,4),D);           % Nnodes x d2
uhd=((wD'*D)\wD')*uD(xx,yy,zz);

%Neumann
x12=x(T.neumann(:,2))-x(T.neumann(:,1)); %x2-x1
y12=y(T.neumann(:,2))-y(T.neumann(:,1)); %y2-y1
z12=z(T.neumann(:,2))-z(T.neumann(:,1)); %z2-z1
x13=x(T.neumann(:,3))-x(T.neumann(:,1)); %x3-x1
y13=y(T.neumann(:,3))-y(T.neumann(:,1)); %y3-y1
z13=z(T.neumann(:,3))-z(T.neumann(:,1)); %z3-z1

Neu_normals=0.5*[y12.*z13-z12.*y13,...
                 z12.*x13-x12.*z13,...
                 x12.*y13-x13.*y12];  % Nneu x 3,  ||n^e||=|e|

xn=formula(:,[1 2 3])*x(T.neumann');  % Nnodes x Nneu
yn=formula(:,[1 2 3])*y(T.neumann');  % Nnodes x Nneu
zn=formula(:,[1 2 3])*z(T.neumann');  % Nnodes x Nneu

qhn=bsxfun(@times,Neu_normals(:,1)',wD'*gx(xn,yn,zn))+...
    bsxfun(@times,Neu_normals(:,2)',wD'*gy(xn,yn,zn))+...
    bsxfun(@times,Neu_normals(:,3)',wD'*gz(xn,yn,zn));  %d2 x Nneu
return
