    function [fh,error]=L2proj3d(f,T,k,formula1,formula2)

%[fh,error]=L2proj3d(f,T,k,formula1,formula2)
%
%Input:
%           f: vectorized function of three variables
%           T: full tetrahedrization
%           k: polynomial degree 
%    formula1: quadrature formula in 3d
%    formula2: quadrature formula in 3d to compute errors
%
%Output:
%          fh: L2 projection of f ; disc P_k function (d3 x Nelts) 
%       error: L^2 error
%
%Last modified: April 2, 2013

P =dubiner3d(2*formula2(:,2)-1,2*formula2(:,3)-1,2*formula2(:,4)-1,k);
wP=bsxfun(@times,formula2(:,5),P);
x=T.coordinates(:,1); x=formula2(:,1:4)*x(T.elements');
y=T.coordinates(:,2); y=formula2(:,1:4)*y(T.elements');
z=T.coordinates(:,3); z=formula2(:,1:4)*z(T.elements');

fh=((wP'*P)\wP')*f(x,y,z);
error=errorElem(T,f,fh,k,formula1);
return