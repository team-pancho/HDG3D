function error=errorFaces(T,p,ph,k,formula)

%error =errorFaces(T,p,ph,k,formula)
%
%Input:
%          T: expanded tetrahedrization
%          p: vectorized function of three variables
%         ph: Pk function on skeleton (d2 x Nfaces)
%          k: polynomial degree
%    formula: quadrature formula in 2d (N x 4 matrix)
%
%Output:
%      error: \| p - ph \|_{h,\partial T_h}
%
%Last modified: March 14, 2013

x=T.coordinates(:,1);x=formula(:,[1 2 3])*x(T.faces(:,[1 2 3])');
y=T.coordinates(:,2);y=formula(:,[1 2 3])*y(T.faces(:,[1 2 3])');
z=T.coordinates(:,3);z=formula(:,[1 2 3])*z(T.faces(:,[1 2 3])');
p=p(x,y,z);  %Nnodes x Nfaces

DB=dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k);
ph=DB*ph;                        %Nnodes x Nfaces

error=sqrt(formula(:,4)'*(p-ph).^2*(T.area).^2);
return