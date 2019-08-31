function Mass=MassMatrix(T,coeffs,k,formula)

%{M1,M2,...}=MassMatrix(T,{c1,c2,...},k,formula)
%Input:
%           T: expanded etrahedrization
% {c1,c2,...}: cell array with vectorized functions of three variables
%           k: polynomial degree
%     formula: quadrature formula in 3d (N x 5 matrix)
%
%Output:
% {M1,M2,...}: each cell is d3 x d3 x Nelts (\int_K c{l} P_i^K P_j^K )           
%
%Last modified: March 14, 2013

Nnodes=size(formula,1);
Nelts=size(T.elements,1);
d3=nchoosek(k+3,3);

x=T.coordinates(:,1);x=formula(:,1:4)*x(T.elements');
y=T.coordinates(:,2);y=formula(:,1:4)*y(T.elements');
z=T.coordinates(:,3);z=formula(:,1:4)*z(T.elements'); % Nnd x Nelts
xhat=formula(:,2);
yhat=formula(:,3);
zhat=formula(:,4);

P=dubiner3d(2*xhat-1,2*yhat-1,2*zhat-1,k);  % Nnd x d3

nMass=size(coeffs,2);
Mass=cell(1,nMass);
for n=1:nMass
    c=coeffs{n};
    C=bsxfun(@times,T.volume',c(x,y,z));           % Nnd x Nelts
    mass=zeros(d3,Nelts*d3);
    for q=1:Nnodes
        mass=mass+kron(C(q,:),formula(q,5)*P(q,:)'*P(q,:));
    end
    mass=reshape(mass,[d3,d3,Nelts]);
    Mass{n}=mass;
end
return