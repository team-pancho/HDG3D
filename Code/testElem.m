function Ints=testElem(f,T,k,formula)

%{Int1,Int2,...}=testElem({f1,f2,...},T,k,formula)
%            Int=testElem({f},T,k,formula)
%
%Input:
% {f1,f2...} : cell array with vectorized functions of three variables
%          T : expanded tetrahedrization
%    formula : 3d quadrature formula (Nnd x 5)
%          k : polynomial degree
%Output:
% {Int1,Int2,...}: each cell is d3 x Nelts (\int_K f{l} P_i^K)
%            Int : matrix \int_K f{l} P_i^K
%Last modified: March 14, 2013

x=T.coordinates(:,1);x=formula(:,1:4)*x(T.elements');
y=T.coordinates(:,2);y=formula(:,1:4)*y(T.elements');
z=T.coordinates(:,3);z=formula(:,1:4)*z(T.elements');
xhat=formula(:,2);
yhat=formula(:,3);
zhat=formula(:,4);

P=dubiner3d(2*xhat-1,2*yhat-1,2*zhat-1,k);
wP=bsxfun(@times,formula(:,5),P);

nInts=size(f,2);
if nInts==1
    Ints=bsxfun(@times,T.volume',wP'*f(x,y,z));
else
    Ints=cell(1,nInts);
    for n=1:nInts
        ff=f{n};
        Ints{n}=bsxfun(@times,T.volume',wP'*ff(x,y,z));
    end
end
return