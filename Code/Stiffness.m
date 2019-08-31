function S=Stiffness(T,k,formula)

% S=stiffness(T,k,formula)
%
% Input:
%   T       : tetrahedrization
%   k       : polynomial degree
%   formula : quadrature formula
% Output:
%   S       : d3 x d3 x Nelts stiffness matrix
% Last modified: March 27, 2013

% Computations in the reference element


xhat=formula(:,2);
yhat=formula(:,3);
zhat=formula(:,4);

[~,px,py,pz]=dubiner3d(2*xhat-1,2*yhat-1,2*zhat-1,k);
px=2*px; wpx=bsxfun(@times,formula(:,5),px);
py=2*py; wpy=bsxfun(@times,formula(:,5),py);
pz=2*pz; wpz=bsxfun(@times,formula(:,5),pz);

Kxx=(1/6)*px'*wpx;
Kxy=(1/6)*px'*wpy;
Kxz=(1/6)*px'*wpz;
Kyy=(1/6)*py'*wpy;
Kyz=(1/6)*py'*wpz;
Kzz=(1/6)*pz'*wpz;

% Geometric constructions

x12=T.coordinates(T.elements(:,2),1)-T.coordinates(T.elements(:,1),1);  
y12=T.coordinates(T.elements(:,2),2)-T.coordinates(T.elements(:,1),2);  
z12=T.coordinates(T.elements(:,2),3)-T.coordinates(T.elements(:,1),3);  
x13=T.coordinates(T.elements(:,3),1)-T.coordinates(T.elements(:,1),1);   
y13=T.coordinates(T.elements(:,3),2)-T.coordinates(T.elements(:,1),2);  
z13=T.coordinates(T.elements(:,3),3)-T.coordinates(T.elements(:,1),3);  
x14=T.coordinates(T.elements(:,4),1)-T.coordinates(T.elements(:,1),1);   
y14=T.coordinates(T.elements(:,4),2)-T.coordinates(T.elements(:,1),2);  
z14=T.coordinates(T.elements(:,4),3)-T.coordinates(T.elements(:,1),3);  
x12=x12'; x13=x13'; x14=x14'; 
y12=y12'; y13=y13'; y14=y14'; 
z12=z12'; z13=z13'; z14=z14'; % row vectors with differences of coords

idetB=1./(x12.*y13.*z14-x12.*y14.*z13-y12.*x13.*z14...
         +y12.*x14.*z13+z12.*x13.*y14-z12.*x14.*y13);

dxdx=kron(idetB.*(y13.*z14-z13.*y14).^2,Kxx)...
     +kron(idetB.*(y13.*z14-z13.*y14).*(-y12.*z14+z12.*y14),Kxy+Kxy')...
     +kron(idetB.*(y13.*z14-z13.*y14).*(y12.*z13-z12.*y13),Kxz+Kxz')...
     +kron(idetB.*(-y12.*z14+z12.*y14).^2,Kyy)...
     +kron(idetB.*(-y12.*z14+z12.*y14).*(y12.*z13-z12.*y13),Kyz+Kyz')...
     +kron(idetB.*(y12.*z13-z12.*y13).^2,Kzz);
dydy=kron(idetB.*(-x13.*z14+z13.*x14).^2,Kxx)...
     +kron(idetB.*(-x13.*z14+z13.*x14).*(x12.*z14-z12.*x14),Kxy+Kxy')...
     +kron(idetB.*(-x13.*z14+z13.*x14).*(-x12.*z13+z12.*x13),Kxz+Kxz')...
     +kron(idetB.*(x12.*z14-z12.*x14).^2,Kyy)...
     +kron(idetB.*(x12.*z14-z12.*x14).*(-x12.*z13+z12.*x13),Kyz+Kyz')...
     +kron(idetB.*(-x12.*z13+z12.*x13).^2,Kzz);
dzdz=kron(idetB.*(x13.*y14-y13.*x14).^2,Kxx)...
     +kron(idetB.*(x13.*y14-y13.*x14).*(-x12.*y14+y12.*x14),Kxy+Kxy')...
     +kron(idetB.*(x13.*y14-y13.*x14).*(x12.*y13-y12.*x13),Kxz+Kxz')...
     +kron(idetB.*(-x12.*y14+y12.*x14).^2,Kyy)...
     +kron(idetB.*(-x12.*y14+y12.*x14).*(x12.*y13-y12.*x13),Kyz+Kyz')...
     +kron(idetB.*(x12.*y13-y12.*x13).^2,Kzz);

 S=dxdx+dydy+dzdz;
 d3=nchoosek(k+3,3);
 Nelts=size(T.elements,1);
 S=reshape(S,[d3 d3 Nelts]);
 
 return
 
 % These matrices are not needed
 
 dxdy=kron(idetB.*(y13.*z14-z13.*y14).*(-x13.*z14+z13.*x14),Kxx)...
     +kron(idetB.*(y13.*z14-z13.*y14).*(x12.*z14-z12.*x14),Kxy)...
     +kron(idetB.*(y13.*z14-z13.*y14).*(-x12.*z13+z12.*x13),Kxz)...
     +kron(idetB.*(-x13.*z14+z13.*x14).*(-y12.*z14+z12.*y14),Kxy')...
     +kron(idetB.*(-y12.*z14+z12.*y14).*(x12.*z14-z12.*x14),Kyy)...
     +kron(idetB.*(-y12.*z14+z12.*y14).*(-x12.*z13+z12.*x13),Kyz)...
     +kron(idetB.*(-x13.*z14+z13.*x14).*(y12.*z13-z12.*y13),Kxz')...
     +kron(idetB.*(x12.*z14-z12.*x14).*(y12.*z13-z12.*y13),Kyz')...
     +kron(idetB.*(y12.*z13-z12.*y13).*(-x12.*z13+z12.*x13),Kzz);
dxdz=kron(idetB.*(y13.*z14-z13.*y14).*(x13.*y14-y13.*x14),Kxx)...
     +kron(idetB.*(y13.*z14-z13.*y14).*(-x12.*y14+y12.*x14),Kxy)...
     +kron(idetB.*(y13.*z14-z13.*y14).*(x12.*y13-y12.*x13),Kxz)...
     +kron(idetB.*(x13.*y14-y13.*x14).*(-y12.*z14+z12.*y14),Kxy')...
     +kron(idetB.*(-y12.*z14+z12.*y14).*(-x12.*y14+y12.*x14),Kyy)...
     +kron(idetB.*(-y12.*z14+z12.*y14).*(x12.*y13-y12.*x13),Kyz)...
     +kron(idetB.*(x13.*y14-y13.*x14).*(y12.*z13-z12.*y13),Kxz')...
     +kron(idetB.*(-x12.*y14+y12.*x14).*(y12.*z13-z12.*y13),Kyz')...
     +kron(idetB.*(y12.*z13-z12.*y13).*(x12.*y13-y12.*x13),Kzz);
dydz=kron(idetB.*(-x13.*z14+z13.*x14).*(x13.*y14-y13.*x14),Kxx)...
     +kron(idetB.*(-x13.*z14+z13.*x14).*(-x12.*y14+y12.*x14),Kxy)...
     +kron(idetB.*(-x13.*z14+z13.*x14).*(x12.*y13-y12.*x13),Kxz)...
     +kron(idetB.*(x13.*y14-y13.*x14).*(x12.*z14-z12.*x14),Kxy')... 
     +kron(idetB.*(x12.*z14-z12.*x14).*(-x12.*y14+y12.*x14),Kyy)...
     +kron(idetB.*(x12.*z14-z12.*x14).*(x12.*y13-y12.*x13),Kyz)...
     +kron(idetB.*(x13.*y14-y13.*x14).*(-x12.*z13+z12.*x13),Kxz')...
     +kron(idetB.*(-x12.*y14+y12.*x14).*(-x12.*z13+z12.*x13),Kyz')...
     +kron(idetB.*(-x12.*z13+z12.*x13).*(x12.*y13-y12.*x13),Kzz);

 