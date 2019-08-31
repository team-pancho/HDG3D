function [fh,error]=L2projskeleton3d(f,T,k,formula1,formula2)

%[fh,error]=L2projskeleton3d(f,T,k,formula1,formula2)
%
%Input:
%           f: vectorized function of three variables
%           T: expanded tetrahedrization
%           k: polynomial degree
%    formula1: quadrature formula in 2d    
%    formula2: quadrature formula in 2d to compute errors
%
%Output:
%          fh: L2 projection on the skeleton (d2 x Nfaces)
%       error: h-weighted error
%
%Last modified: April 2, 2013

D =dubiner2d(2*formula1(:,2)-1,2*formula1(:,3)-1,k);   %Nnodes x d2
wD=bsxfun(@times,formula1(:,4),D);
x=T.coordinates(:,1); x=formula1(:,1:3)*x(T.faces(:,1:3)');
y=T.coordinates(:,2); y=formula1(:,1:3)*y(T.faces(:,1:3)');
z=T.coordinates(:,3); z=formula1(:,1:3)*z(T.faces(:,1:3)');

fh=((wD'*D)\wD')*f(x,y,z);

error=errorFaces(T,f,fh,k,formula2);
return