function [C,Cf,A1,A2,Af]=localsolvers3d(km,c,f,tau,T,k,formulas)

%[C,Cf,A1,A2,Af]=localsolvers3d(km,c,f,tau,T,k,{for1,for2,for3})
%
%Input: 
%      km, c, f: vectorized functions of three variables 
%           tau: penalization parameter for HDG (Nelts x 4)
%             T: expanded tetrahedrization
%             k: polynomial degree
%          for1: quadrature formula 3d (for mass matrix)
%          for2: quadrature formula 3d (for convection matrices)
%          for3: quadrature formula 2d
%
%Output:
%             C:   4*d2 x 4*d2 x Nelts
%            Cf:   4*d2 x Nelts
%            A1:   4*d3 x 4*d3 x Nelts
%            A2:   4*d3 x 4*d2 x Nelts
%            Af:   4*d3 x Nelts
%
%Last modified: Nov 29, 2012

Nelts=size(T.elements,1);
d2=nchoosek(k+2,2);    
d3=nchoosek(k+3,3); 

f=testElem(f,T,k,formulas{1});
Af=zeros(4*d3,Nelts);
Af(3*d3+1:4*d3,:)=f;

Mass=MassMatrix(T,{km,c},k,formulas{1}); 
Mk=Mass{1};Mc=Mass{2};
[Cx,Cy,Cz]=ConvMatrix(T,k,formulas{2});
[tauPP,tauDP,nxDP,nyDP,nzDP,tauDD]=matricesFace(T,tau,k,formulas{3});

A1=zeros(4*d3,4*d3,Nelts);
A2=zeros(4*d3,4*d2,Nelts);
C=zeros(4*d2,4*d2,Nelts);
Cf=zeros(4*d2,Nelts);

O=zeros(d3,d3,Nelts);
A1=[Mk        ,O        ,O        ,-permute(Cx,[2 1 3]);...
    O         ,Mk       ,O        ,-permute(Cy,[2 1 3]);...
    O         ,O        ,Mk       ,-permute(Cz,[2 1 3]);...
    Cx        ,Cy       ,Cz       ,Mc+tauPP];
A2=[permute(nxDP,[2 1 3]);...
    permute(nyDP,[2 1 3]);...
    permute(nzDP,[2 1 3]);...
    -permute(tauDP,[2 1 3])];            

parfor i=1:Nelts            
    C(:,:,i)=[nxDP(:,:,i) nyDP(:,:,i) nzDP(:,:,i) tauDP(:,:,i)]/A1(:,:,i)*A2(:,:,i)...
              +tauDD(:,:,i);
    Cf(:,i) =[nxDP(:,:,i) nyDP(:,:,i) nzDP(:,:,i) tauDP(:,:,i)]/A1(:,:,i)*Af(:,i);
end
return