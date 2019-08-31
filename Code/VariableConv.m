function [convx,convy,convz]=VariableConv(T,vx,vy,vz,k,formula)

%[convx,convy,convz]=VariableConv(T,vx,vy,vz,k,formula)
%
%Input:
%         T: expanded tetrahedrization
%  vx,vy,vz: vectorized functions of three variables
%         k: polynomial degree
%   formula: quadrature formula in 3d (N x 5 matrix)
%
%Output:
%     convx: d3 x d3 x Nelts ( \int_K vx P_i^K \partial_x P_j^K )
%     convy: d3 x d3 x Nelts ( \int_K vy P_i^K \partial_y P_j^K )
%     convz: d3 x d3 x Nelts ( \int_K vz P_i^K \partial_z P_j^K )
%
%Last modified: May 31, 2012

Nelts=size(T.elements,1);
Nnodes=size(formula,1);
d3=nchoosek(k+3,3);

x=T.coordinates(:,1);x=x(T.elements)*formula(:,1:4)'; %Nelts x Nnodes
y=T.coordinates(:,2);y=y(T.elements)*formula(:,1:4)';
z=T.coordinates(:,3);z=z(T.elements)*formula(:,1:4)'; 

convx=zeros(d3,Nelts*d3);
convy=zeros(d3,Nelts*d3);
convz=zeros(d3,Nelts*d3);

xhat=formula(:,2);
yhat=formula(:,3);
zhat=formula(:,4);

[p,px,py,pz]=dubiner3d(2*xhat-1,2*yhat-1,2*zhat-1,k); % Nnodes x d3
px=2*px;
py=2*py;
pz=2*pz;

x12=T.coordinates(T.elements(:,2),1)-T.coordinates(T.elements(:,1),1); %x2-x1
x13=T.coordinates(T.elements(:,3),1)-T.coordinates(T.elements(:,1),1); %x3-x1
x14=T.coordinates(T.elements(:,4),1)-T.coordinates(T.elements(:,1),1); %x4-x1
y12=T.coordinates(T.elements(:,2),2)-T.coordinates(T.elements(:,1),2); %y2-y1
y13=T.coordinates(T.elements(:,3),2)-T.coordinates(T.elements(:,1),2); %y3-y1
y14=T.coordinates(T.elements(:,4),2)-T.coordinates(T.elements(:,1),2); %y4-y1
z12=T.coordinates(T.elements(:,2),3)-T.coordinates(T.elements(:,1),3); %z2-z1
z13=T.coordinates(T.elements(:,3),3)-T.coordinates(T.elements(:,1),3); %z3-z1
z14=T.coordinates(T.elements(:,4),3)-T.coordinates(T.elements(:,1),3); %z4-z1

a=y13.*z14-y14.*z13;
b=y14.*z12-y12.*z14;
c=y12.*z13-y13.*z12;
d=x14.*z13-x13.*z14;
e=x12.*z14-x14.*z12;
f=x13.*z12-x12.*z13;
g=x13.*y14-x14.*y13;
h=x14.*y12-x12.*y14;
i=x12.*y13-x13.*y12;

A=bsxfun(@times,a,vx(x,y,z));
B=bsxfun(@times,b,vx(x,y,z));
C=bsxfun(@times,c,vx(x,y,z));
D=bsxfun(@times,d,vy(x,y,z));
E=bsxfun(@times,e,vy(x,y,z));
F=bsxfun(@times,f,vy(x,y,z));
G=bsxfun(@times,g,vz(x,y,z));
H=bsxfun(@times,h,vz(x,y,z));
I=bsxfun(@times,i,vz(x,y,z));

for i=1:Nnodes
    convx=convx+kron(A(:,i)',1/6*formula(i,5)*p(i,:)'*px(i,:))+...
                kron(B(:,i)',1/6*formula(i,5)*p(i,:)'*py(i,:))+...
                kron(C(:,i)',1/6*formula(i,5)*p(i,:)'*pz(i,:));
    convy=convy+kron(D(:,i)',1/6*formula(i,5)*p(i,:)'*px(i,:))+...
                kron(E(:,i)',1/6*formula(i,5)*p(i,:)'*py(i,:))+...
                kron(F(:,i)',1/6*formula(i,5)*p(i,:)'*pz(i,:));
    convz=convz+kron(G(:,i)',1/6*formula(i,5)*p(i,:)'*px(i,:))+...
                kron(H(:,i)',1/6*formula(i,5)*p(i,:)'*py(i,:))+...
                kron(I(:,i)',1/6*formula(i,5)*p(i,:)'*pz(i,:));
    
end
convx=reshape(convx,[d3,d3,Nelts]);
convy=reshape(convy,[d3,d3,Nelts]);
convz=reshape(convz,[d3,d3,Nelts]);
return