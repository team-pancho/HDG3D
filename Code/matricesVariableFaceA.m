function mat=matricesVariableFaceA(T,al,pw,k,formula)

% {m1,m2,...}=matricesFaceVariable(T,{a1,a2,...},{pw1,pw2...},k,form)
%
% Input:
%      T           : full tetrahedrization
%      a1,a2,...   : vectorized function of three variables
%      pw1,pw2,... : piecewise constant function on skeleton 4 x Nfaces
%      k           : polynomial degree
%      form        : 2-dimensional quadrature rule
% Output:
%      m1,m2,...   : 4d2 x d3 x Nelts
%
% Last modified: April 5, 2013

% Evaluations of Dubiner3d functions on all faces

Nelts=size(T.elements,1);
Nnodes=size(formula,1);
T.perm=T.perm';                       % 4 x Nelts

d2=nchoosek(k+2,2);   
d3=nchoosek(k+3,3);   

s=formula(:,2);
t=formula(:,3);
O=zeros(size(s));
points3d=[s,t,O;...
          s,O,t;...
          O,s,t;...
          s,t,1-s-t];
                 
P =dubiner3d(2*points3d(:,1)-1,2*points3d(:,2)-1,2*points3d(:,3)-1,k);
wP=bsxfun(@times,[formula(:,4);...
                  formula(:,4);...
                  formula(:,4);...
                  formula(:,4)],P);
              
% Evaluations of alpha on all faces

x=T.coordinates(:,1);
y=T.coordinates(:,2);
z=T.coordinates(:,3);
points3d=[1-sum(points3d,2) points3d];  % barycentric coordinates
x=points3d*x(T.elements');
y=points3d*y(T.elements');
z=points3d*z(T.elements');

% Evaluations of Dubiner2d polynomials on six different configurations

points2d=[s,t;...
          t,s;...
          1-s-t,s;...
          s,1-s-t;...
          t,1-s-t;...
          1-s-t,t];
D=dubiner2d(2*points2d(:,1)-1,2*points2d(:,2)-1,k);    % 6*Nnodes x d2

% Construction of the matrix

for c=1:length(al)
    alpha=al{c};
    alpha=alpha(x,y,z);       % 4Nnodes X Nelts
    pwct =pw{c};
    matrix=zeros(4*d2,d3*Nelts);
    for l=1:4
        rows=1+(l-1)*d2:(l*d2);
        for mu=1:6
            pwctmu=pwct(l,:).*(T.perm(l,:)==mu);
            for r=1:Nnodes
                rowP=(l-1)*Nnodes+r;
                rowD=(mu-1)*Nnodes+r;
                matrix(rows,:)=matrix(rows,:)...
                    +kron(pwctmu.*alpha(rowP,:),D(rowD,:)'*wP(rowP,:));
            end
        end
    end
    mat{c}=reshape(matrix,[4*d2,d3,Nelts]);
end

              
              

