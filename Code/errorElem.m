function error=errorElem(T,p,ph,k,formula)

%error=errorElem(T,p,ph,k,formula)
%
%Input:
%          T: expanded tetrahedrization
%          p: vectorized function of three variables
%         ph: discontinuous Pk function (d3 x Nelts)
%          k: polynomial degree
%    formula: quadrature formula in 3d (N x 5 matrix)
%
%Output:
%      error: \| p - ph \|_{L^2}
%
%Last modified: May 31, 2012

x=T.coordinates(:,1); x=formula(:,1:4)*x(T.elements');
y=T.coordinates(:,2); y=formula(:,1:4)*y(T.elements');
z=T.coordinates(:,3); z=formula(:,1:4)*z(T.elements');
p=p(x,y,z);

xhat=formula(:,2);
yhat=formula(:,3);
zhat=formula(:,4);
B=dubiner3d(2*xhat-1,2*yhat-1,2*zhat-1,k);
ph=B*ph;

error=sqrt(formula(:,5)'*(p-ph).^2*T.volume);
return