function mat=matricesVariableFaceB(T,al,pw,k,formula)

% {m1,m2,...}=matricesFaceVariable(T,{a1,a2,...},{pw1,pw2,...},k,form)
%
% Input:
%      T           : full tetrahedrization
%      a1,a2,...   : vectorized function of three variables
%      pw1,pw2,... : piecewise constant function on skeleton 4 x Nfaces
%      k           : polynomial degree
%      form        : 2-dimensional quadrature rule
% Output:
%      m1,m2,...   : 4d2 x 4d2 x Nelts
%
% Last modified: April 2, 2013

% Parameters

Nelts=size(T.elements,1);
Nnodes=size(formula,1);
d2=nchoosek(k+2,2);   
O=zeros(d2,d2,Nelts);
              
% Evaluations of alpha and Dubiner 2d on all faces

x=T.coordinates(:,1);x=formula(:,1:3)*x(T.faces(:,1:3)');
y=T.coordinates(:,2);y=formula(:,1:3)*y(T.faces(:,1:3)');
z=T.coordinates(:,1);z=formula(:,1:3)*z(T.faces(:,1:3)');
D=dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k);

% Construction of the matrix

for c=1:length(al);
    alpha=al{c};
    alpha=alpha(x,y,z);       % Nnodes X Nfaces
    pwct=pw{c};
    for l=1:4
        dgm{l}=zeros(d2,d2*Nelts);
        for r=1:Nnodes
            dgm{l}=dgm{l}+...
                kron(pwct(l,:).*alpha(r,T.facebyele(:,l)),...
                     formula(r,4)*D(r,:)'*D(r,:));
        end
        dgm{l}=reshape(dgm{l},[d2,d2,Nelts]);
    end
    mat{c}=[dgm{1} O      O      O      ;...
            O      dgm{2} O      O      ;...
            O      O      dgm{3} O      ;...
            O      O      O      dgm{4}];
end
return

    
    
              
              

