function [tauPP,tauDP,nxDP,nyDP,nzDP,tauDD]=matricesFace(T,tau,k,formula)

%[tauPP,tauDP,nxDP,nyDP,nzDP,tauDD]=matricesFace(T,tau,k,formula)
%
%Input:
%             T: expanded terahedrization
%           tau: penalization parameter for HDG (Nelts x 4)
%             k: polynomial degree 
%       formula: quadrature formula in 2d (N x 4 matrix)
%Output:
%        tauPP :   d3 x d3   x Nelts, with <tau P_i,P_j>_{\partial K}
%        tauDP : 4*d2 x d3   x Nelts, with <tau D_i,P_j>_e, e\in E(K)
%         nxDP : 4*d2 x d3   x Nelts, with <nx D_i,P_j>_e, e\in E(K)
%         nyDP : 4*d2 x d3   x Nelts, with <ny D_i,P_j>_e, e\in E(K)
%         nzDP : 4*d2 x d3   x Nelts, with <nz D_i,P_j>_e, e\in E(K)
%        tauDD : 4*d2 x 4*d2 x Nelts, block diag <tau D_i, D_j>_e, e\in E(K)
%
%Last modified: March 14, 2013

Nelts=size(T.elements,1);
Nnodes=size(formula,1);
TauArea=T.area(T.facebyele').*tau;    % 4 x Nelts
T.perm=T.perm';                       % 4 x Nelts

d2=nchoosek(k+2,2);   
d3=nchoosek(k+3,3);   

s=formula(:,2);
t=formula(:,3);
weights=formula(:,4);

% Computation <tau*Pi,Pj>

O=zeros(size(s));
points3d=[s,t,O;...
          s,O,t;...
          O,s,t;...
          s,t,1-s-t];
                 
pb=dubiner3d(2*points3d(:,1)-1,2*points3d(:,2)-1,2*points3d(:,3)-1,k);  %4*Nnodes x d3
pbweights=bsxfun(@times,[formula(:,4);...
                         formula(:,4);...
                         formula(:,4);...
                         formula(:,4)],pb);
pbpb1=pbweights(1:Nnodes,:)'*pb(1:Nnodes,:);
pbpb2=pbweights(Nnodes+1:2*Nnodes,:)'*pb(Nnodes+1:2*Nnodes,:);
pbpb3=pbweights(2*Nnodes+1:3*Nnodes,:)'*pb(2*Nnodes+1:3*Nnodes,:);
pbpb4=pbweights(3*Nnodes+1:4*Nnodes,:)'*pb(3*Nnodes+1:4*Nnodes,:);

tauPP=kron(TauArea(1,:),pbpb1)+kron(TauArea(2,:),pbpb2)...
      +kron(TauArea(3,:),pbpb3)+kron(TauArea(4,:),pbpb4);
tauPP=reshape(tauPP,[d3,d3,Nelts]);

% Computation <alpha*D,P>, alpha=tau,nx,ny,nz, 

pb=[pb(1:Nnodes,:),pb(Nnodes+1:2*Nnodes,:),...
    pb(2*Nnodes+1:3*Nnodes,:),pb(3*Nnodes+1:4*Nnodes,:)];  % Nnodes x 4*d3
points2d=[s,t;...
          t,s;...
          1-s-t,s;...
          s,1-s-t;...
          t,1-s-t;...
          1-s-t,t];
db=dubiner2d(2*points2d(:,1)-1,2*points2d(:,2)-1,k);    % 6*Nnodes x d2
db=[db(1:Nnodes,:),db(Nnodes+1:2*Nnodes,:),...
    db(2*Nnodes+1:3*Nnodes,:),db(3*Nnodes+1:4*Nnodes,:),...
    db(4*Nnodes+1:5*Nnodes,:),db(5*Nnodes+1:6*Nnodes,:)]; % Nnodes x 6*d2
db=bsxfun(@times,weights,db);
allproducts=db'*pb;             %6*d2 x 4*d3

block2=@(x) (1+(x-1)*d2):(x*d2);
block3=@(x) (1+(x-1)*d3):(x*d3);

tauDP=zeros(4*d2,d3*Nelts);
nxDP=zeros(4*d2,d3*Nelts);
nyDP=zeros(4*d2,d3*Nelts);
nzDP=zeros(4*d2,d3*Nelts);
for l=1:4
    Nx=T.normals(:,3*(l-1)+1)';
    Ny=T.normals(:,3*(l-1)+2)';
    Nz=T.normals(:,3*(l-1)+3)';
    for mu=1:6
        taumu=TauArea(l,:).*(T.perm(l,:)==mu);
        tauDP(block2(l),:)=tauDP(block2(l),:)+...
            kron(taumu,allproducts(block2(mu),block3(l)));
        
        nxmu=Nx.*(T.perm(l,:)==mu);
        nxDP(block2(l),:)=nxDP(block2(l),:)+...
            kron(nxmu,allproducts(block2(mu),block3(l)));
        
        nymu=Ny.*(T.perm(l,:)==mu);
        nyDP(block2(l),:)=nyDP(block2(l),:)+...
            kron(nymu,allproducts(block2(mu),block3(l)));
        
        nzmu=Nz.*(T.perm(l,:)==mu);
        nzDP(block2(l),:)=nzDP(block2(l),:)+...
            kron(nzmu,allproducts(block2(mu),block3(l)));
    end
end
tauDP=reshape(tauDP,[4*d2,d3,Nelts]);
nxDP=reshape(nxDP,[4*d2,d3,Nelts]);
nyDP=reshape(nyDP,[4*d2,d3,Nelts]);
nzDP=reshape(nzDP,[4*d2,d3,Nelts]);

% Computation tauDD 

d=dubiner2d(2*s-1,2*t-1,k);
dweights=bsxfun(@times,d,weights);
dwd=dweights'*d;
tauDD=zeros(4*d2,4*d2,Nelts);
for l=1:4
    tauDD(block2(l),block2(l),:)=reshape(kron(TauArea(l,:),dwd),...
                                         [d2,d2,Nelts]);
end
return