function [Pqx,Pqy,Pqz,Pu]=projectHDG3d(T,coeffs,k,tau,formulas)

%[Pqx,Pqy,Pqz,Pu]=projectHDG3d(T,coeffs,k,tau,formulas)
%
%Input:
%          T: Tetrahedrization
%     coeffs: each cell is a vectorized function of three variables
%          k: degree of polynomial
%   formulas: quadrature matrix for mass matrix and convection matrix
%
%Output:
%        Pqx: \Pi qx, 
%        Pqy: \Pi qy, 
%        Pqz: \Pi qz, 
%        Pu: \Pi u, 

%Last modified: April 9, 2013

d2=nchoosek(k+2,2); 
d3=nchoosek(k+3,3); %Nbasis
if k>0
    d4=nchoosek(k+2,3); %dim P_{k-1}(K)
else
    d4=0;
end 
block3=@(x) (1+(x-1)*d3):(x*d3);
a=@(x,y,z) 1+0.*x;

Mass=MassMatrix(T,{a},k,formulas{1}); 

Mass=Mass{1};
rows=1+d4:d3;
Mass(rows,:,:)=[];

[tauPP,tauDP,nxDP,nyDP,nzDP]=matricesFace(T,tau,k,formulas{3});
Nelts=size(T.elements,1);
O=zeros(d4,d3,Nelts);
M=[Mass O    O    O   ;
   O    Mass O    O   ;
   O    O    Mass O   ;
   O    O    O    Mass; 
   nxDP nyDP nzDP tauDP];

if k>0
    Ints=testElem(coeffs,T,k-1,formulas{2});
    qxP=Ints{1};
    qyP=Ints{2}; 
    qzP=Ints{3}; 
    uP=Ints{4};
else
    qxP=zeros(0,Nelts);
    qyP=zeros(0,Nelts);
    qzP=zeros(0,Nelts);
    uP=zeros(0,Nelts);
end

Ints=testFaces(coeffs,T,k,formulas{3});
qxD=Ints{1}; 
qyD=Ints{2}; 
qzD=Ints{3}; 
 ud=Ints{4};

nx=T.normals(:,[1 4 7 10]);
ny=T.normals(:,[2 5 8 11]);
nz=T.normals(:,[3 6 9 12]);
nx=nx./T.area(T.facebyele);
ny=ny./T.area(T.facebyele);
nz=nz./T.area(T.facebyele);

QxD=zeros(d2,4,Nelts);
QyD=QxD;
QzD=QxD;
UD =QxD;

Pqx=zeros(d3,Nelts);
Pqy=Pqx;
Pqz=Pqx;
Pu =Pqx;

parfor i=1:Nelts
    QxD(:,:,i)=bsxfun(@times,qxD(:,T.facebyele(i,:)),nx(i,:));
    QyD(:,:,i)=bsxfun(@times,qyD(:,T.facebyele(i,:)),ny(i,:));
    QzD(:,:,i)=bsxfun(@times,qzD(:,T.facebyele(i,:)),nz(i,:));
     UD(:,:,i)=bsxfun(@times,ud(:,T.facebyele(i,:)),tau(:,i)');
    qxv=QxD(:,:,i);
    qyv=QyD(:,:,i);
    qzv=QzD(:,:,i);
     uv=UD(:,:,i);
    
    rhs=[qxP(:,i);
         qyP(:,i);
         qzP(:,i);
          uP(:,i);
          qxv(:)+qyv(:)+qzv(:)+uv(:)];
     vect=M(:,:,i)\rhs;
      
     Pqx(:,i)=vect(block3(1));
     Pqy(:,i)=vect(block3(2));
     Pqz(:,i)=vect(block3(3));
      Pu(:,i)=vect(block3(4));
end
return