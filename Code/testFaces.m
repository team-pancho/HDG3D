function Ints=testFaces(f,T,k,formula)

%{Int1,Int2,...}=testFaces({f1,f2,...},T,k,formula)
%            Int=testFaces({f1},T,k,formula)
%
%Input:
%  {f1,f2,...}: each cell is a vectorized function of three variables
%            T: expanded tetrahedrization
%      formula: quadrature formula in 2d (N x 4 matrix)
%            k: polynomial
%
%Output:
%  {Int1,Ints2,...}: each cell is a d2 x Nfaces matrix
%                    ( \int_e f{l} D_i^e )
%
%Last modified: March 14, 2013

x=T.coordinates(:,1); x=formula(:,1:3)*x(T.faces(:,[1 2 3])');
y=T.coordinates(:,2); y=formula(:,1:3)*y(T.faces(:,[1 2 3])');
z=T.coordinates(:,3); z=formula(:,1:3)*z(T.faces(:,[1 2 3])');

DB=dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k);
DB=bsxfun(@times,formula(:,4),DB);

nInts=size(f,2);
if nInts==1
    Ints=bsxfun(@times,T.area',DB'*f(x,y,z));
else
    Ints=cell(1,nInts);
    for n=1:nInts
        ff=f{n};
        Ints{n}=bsxfun(@times,T.area',DB'*ff(x,y,z));
    end
end
return