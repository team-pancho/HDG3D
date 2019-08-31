function [convx,convy,convz]=ConvMatrix(T,k,formula)

%[convx,convy,convz]=ConvMatrix(T,k,formula)
%
%Input:
%        T: expanded tetrahedrization
%        k: polynomial degree
%  formula: quadrature formula in 3d (N x 5 matrix)
%
%Output:
%    convx: d3 x d3 x Nelts  ( \int_K P_i^K \partial_x P_j^K )
%    convy: d3 x d3 x Nelts  ( \int_K P_i^K \partial_y P_j^K )
%    convz: d3 x d3 x Nelts  ( \int_K P_i^K \partial_z P_j^K )
%
%Last modified: March 14, 2013

Nelts=size(T.elements,1);
d3=nchoosek(k+3,3);

xhat=formula(:,2);
yhat=formula(:,3);
zhat=formula(:,4);

[p,px,py,pz]=dubiner3d(2*xhat-1,2*yhat-1,2*zhat-1,k);
px=2*px;
py=2*py;
pz=2*pz;
wp=bsxfun(@times,formula(:,5),p);

convhatx=1/6*wp'*px;
convhaty=1/6*wp'*py;
convhatz=1/6*wp'*pz;

x12=T.coordinates(T.elements(:,2),1)-T.coordinates(T.elements(:,1),1); %x2-x1
x13=T.coordinates(T.elements(:,3),1)-T.coordinates(T.elements(:,1),1); %x3-x1
x14=T.coordinates(T.elements(:,4),1)-T.coordinates(T.elements(:,1),1); %x4-x1
y12=T.coordinates(T.elements(:,2),2)-T.coordinates(T.elements(:,1),2); %y2-y1
y13=T.coordinates(T.elements(:,3),2)-T.coordinates(T.elements(:,1),2); %y3-y1
y14=T.coordinates(T.elements(:,4),2)-T.coordinates(T.elements(:,1),2); %y4-y1
z12=T.coordinates(T.elements(:,2),3)-T.coordinates(T.elements(:,1),3); %z2-z1
z13=T.coordinates(T.elements(:,3),3)-T.coordinates(T.elements(:,1),3); %z3-z1
z14=T.coordinates(T.elements(:,4),3)-T.coordinates(T.elements(:,1),3); %z4-z1

axx=y13.*z14-y14.*z13;
axy=y14.*z12-y12.*z14;
axz=y12.*z13-y13.*z12;
ayx=x14.*z13-x13.*z14;
ayy=x12.*z14-x14.*z12;
ayz=x13.*z12-x12.*z13;
azx=x13.*y14-x14.*y13;
azy=x14.*y12-x12.*y14;
azz=x12.*y13-x13.*y12;

convx=kron(axx',convhatx)+kron(axy',convhaty)+kron(axz',convhatz);
convy=kron(ayx',convhatx)+kron(ayy',convhaty)+kron(ayz',convhatz);
convz=kron(azx',convhatx)+kron(azy',convhaty)+kron(azz',convhatz);

convx=reshape(convx,[d3,d3,Nelts]);
convy=reshape(convy,[d3,d3,Nelts]);
convz=reshape(convz,[d3,d3,Nelts]);
return