function [Uh,Qxh,Qyh,Qzh,Uhat,system,solvers]=...
          HDG3dCD(km,c,f,beta,tau,T,k,formulas,uD,gx,gy,gz,varargin)

%[Uh,Qxh,Qyh,Qzh,Uhat]=HDG3dCD(km,c,f,{bx,by,bz},...
%                       tau,T,k,formulas,uD,gx,gy,gz)
%[Uh,Qxh,Qyh,Qzh,Uhat]=HDG3dCD(km,c,f,{bx,by,bz},...
%                       tau,T,k,formulas,uD,gx,gy,gz,0)
%[~,~,~,~,~,system,solvers]=HDG3dCD(km,c,f,{bx,by,bz},...
%                       tau,T,k,formulas,uD,gx,gy,gz,1)
%
%Input:
%        km   : vectorized function (kappa^{-1}; kappa=diffusion parameter)
%        c    : vectorized function (reaction parameter)
%        f    : vectorized function (source term)
%  {bx,by,bz} : vectorized functions (convection field)
%        tau  : penalization parameter for HDG (Nelts x 4)
%        T    : expanded tetrahedrization
%        k    : polynomial degree
%     formulas: {for1,for2,for3,for4} 
%               (quadrature formulas as used by localsolvers3dCD)
%        uD   : Dirichlet data; vectorized function
%    gx,gy,gz : Neumann data (corresponds to kappa*grad(u)); vectorized fns  
%
%Output:
%      Uh     : d3 x Nelts,   matrix with uh
%      Qxh    : d3 x Nelts,   matrix with qxh
%      Qyh    : d3 x Nelts,   matrix with qyh
%      Qzh    : d3 x Nelts,   matrix with qzh
%      Uhat   : d2 x Nelts    matrix with uhhat
%    system   : {HDGmatrix,HDGrhs,list of free d.o.f.,list of dir d.o.f.}
%    solvers  : {A1,A2,Af}    local solvers
%
%Last modified: April 11, 2013

if nargin==13
    export=varargin{1};
else
    export=0;
end

d2=nchoosek(k+2,2);    
d3=nchoosek(k+3,3); 
block3=@(x) (1+(x-1)*d3):(x*d3);
Nelts =size(T.elements,1);
Nfaces=size(T.faces,1);
Ndir  =size(T.dirichlet,1);
Nneu  =size(T.neumann,1);

%Matrices for assembly process

face=T.facebyele';    % 4 x Nelts
face=(face(:)-1)*d2;  % First degree of freedom of each face by element                                               
face=bsxfun(@plus,face,1:d2);    %4*Nelts x d2 (d.o.f. for each face)                                              
face=reshape(face',4*d2,Nelts);  %d.o.f. for the 4 faces of each element

[J,I]=meshgrid(1:4*d2);
R=face(I(:),:); R=reshape(R,4*d2,4*d2,Nelts);
C=face(J(:),:); C=reshape(C,4*d2,4*d2,Nelts); 
      % R_ij^K d.o.f. for local (i,j) d.o.f. in element K ; R_ij^K=C_ji^K
RowsRHS=reshape(face,4*d2*Nelts,1);

dirfaces=(T.dirfaces(:)-1)*d2;
dirfaces=bsxfun(@plus,dirfaces,1:d2);
dirfaces=reshape(dirfaces',d2*Ndir,1);

free=((1:Nfaces)'-1)*d2;
free=bsxfun(@plus,free,1:d2);
free=reshape(free',d2*Nfaces,1);
free(dirfaces)=[];

neufaces=(T.neufaces(:)-1)*d2;
neufaces=bsxfun(@plus,neufaces,1:d2);
neufaces=reshape(neufaces',d2*Nneu,1);

%Local solvers and global system

[M1,Cf,A1,A2,Af]=localsolvers3dCD(km,c,f,beta,tau,T,k,formulas);
M=sparse(R(:),C(:),M1(:));
phif=accumarray(RowsRHS,Cf(:));

[uhatD,qhatN]=BC3d(uD,gx,gy,gz,T,k,formulas{3});

%Dirichlet BC
Uhatv=zeros(d2*Nfaces,1);
Uhatv(dirfaces)=uhatD;       %uhat stored as a vector: d2*Nfaces

%RHS
rhs=zeros(d2*Nfaces,1);
rhs(free)=phif(free);
qhatN=reshape(qhatN,d2*Nneu,1);   % qhatN stored as a vector: d2*Nneu
rhs(neufaces)=rhs(neufaces)+qhatN;
rhs=rhs-M(:,dirfaces)*Uhatv(dirfaces);

if export
    system={M,rhs,free,dirfaces};
    solvers={A1,A2,Af};
    Uh=[];
    Qxh=[];
    Qyh=[];
    Qzh=[];
    Uhat=[];
    return
else
    system=[];
    solvers=[];
end

Uhatv(free)=M(free,free)\rhs(free);
Uhat=reshape(Uhatv,d2,Nfaces);

%Uh Qxh Qyh Qzh

faces=T.facebyele'; faces=faces(:);
uhhataux=reshape(Uhat(:,faces),[4*d2,Nelts]);
sol=zeros(4*d3,Nelts);
parfor K=1:Nelts
    sol(:,K)=A1(:,:,K)\(Af(:,K)-A2(:,:,K)*uhhataux(:,K));
end

Qxh=sol(block3(1),:);
Qyh=sol(block3(2),:);
Qzh=sol(block3(3),:);
Uh =sol(block3(4),:);
return