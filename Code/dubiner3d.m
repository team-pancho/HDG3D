function [db,dbx,dby,dbz]=dubiner3d(x,y,z,k)
% function [db,dbx,dby,dbz]=dubiner3d(x,y,z,k)
% input:
%     x,y,z: the coordinates at which we evaluate the Dubiner basis
%         k: the degree of the polynomial. Must be a column vector.
% output:
%       db: the value of dubiner basis at given coordinates, Nnodes x Nbasis 
%      dbx: the derivative to x of the Dubiner basis, Nnodes x Nbasis  
%      dby: the derivative to y of the Dubiner basis, Nnodes x Nbasis
%      dbz: the derivative to y of the Dubiner basis, Nnodes x Nbasis
% Last modified: May 16 2012

Nnodes=size(x,1);
Nbasis=nchoosek(k+3,3);
db=zeros(Nnodes,Nbasis);
dbx=zeros(Nnodes,Nbasis);
dby=zeros(Nnodes,Nbasis);
dbz=zeros(Nnodes,Nbasis);
a=@(al,be,n) ( (2*n+1+al+be)*(2*n+2+al+be) )/ ( 2*(n+1)*(n+1+al+be) );
b=@(al,be,n) ( (al^2-be^2)*(2*n+1+al+be) ) / ( 2*(n+1)*(2*n+al+be)*...
                                                          (n+1+al+be)  );
c=@(al,be,n) ( (n+al)*(n+be)*(2*n+2+al+be) ) / ( (n+1)*(n+1+al+be)*...
                                                (2*n+al+be)  );

% Locate the basis with index p,q,r
index=zeros(k+1,k+1,k+1);
count=1;
%
for j=0:k
    for i1=j:-1:0
        for i2=j-i1:-1:0
            index(i1+1,i2+1,j-i1-i2+1)=count;
            count=count+1;
        end
    end
end
loc=@(i1,i2,i3)  index(i1+1,i2+1,i3+1);


F1= (2+2*x+y+z)/2;  F1x= 1+0*x;   F1y= 1/2+0*x;    F1z= 1/2+0*x;
F2= ((y+z)/2).^2;   F2x= 0.*x;    F2y= 1/2*(y+z);  F2z= 1/2*(y+z); 
F3= (2+3*y+z)/2;    F3x= 0+0.*x;  F3y= 3/2+0.*x;   F3z= 1/2+0.*x;
F4= (1+2*y+z)/2;    F4x= 0+0.*x;  F4y= 1+0.*x;     F4z= 1/2+0.*x;
F5= (1-z)/2;        F5x= 0+0.*x;  F5y= 0+0.*x;     F5z= -1/2+0.*x;

db(:,loc(0,0,0))=1;
dbx(:,loc(0,0,0))=0*x;   dby(:,loc(0,0,0))=0*y;  dbz(:,loc(0,0,0))=0*z;
if k>0
    db(:,loc(1,0,0))=F1;
    dbx(:,loc(1,0,0))=F1x; dby(:,loc(1,0,0))=F1y; dbz(:,loc(1,0,0))=F1z;
end
for p=1:k-1
    db(:,loc(p+1,0,0))=(2*p+1)/(p+1)*F1.*db(:,loc(p,0,0))...
        -(p/(p+1))*F2.*db(:,loc(p-1,0,0));
    dbx(:,loc(p+1,0,0))=(2*p+1)/(p+1)*( F1x.*db(:,loc(p,0,0)) + ...
        F1.*dbx(:,loc(p,0,0)) ) - (p/(p+1))*(  F2x.*db(:,loc(p-1,0,0))+ ...
        F2.*dbx(:,loc(p-1,0,0))  );
    dby(:,loc(p+1,0,0))=(2*p+1)/(p+1)*( F1y.*db(:,loc(p,0,0)) + ...
        F1.*dby(:,loc(p,0,0)) ) - (p/(p+1))*(  F2y.*db(:,loc(p-1,0,0))+ ...
        F2.*dby(:,loc(p-1,0,0))  );
    dbz(:,loc(p+1,0,0))=(2*p+1)/(p+1)*( F1z.*db(:,loc(p,0,0)) + ...
        F1.*dbz(:,loc(p,0,0)) ) - (p/(p+1))*(  F2z.*db(:,loc(p-1,0,0))+ ...
        F2.*dbz(:,loc(p-1,0,0))  );
end

for p=0:k-1
    db(:,loc(p,1,0))=(p*(1+y)+F3).*db(:,loc(p,0,0));
    dbx(:,loc(p,1,0)) = F3x.*db(:,loc(p,0,0)) +...
        dbx(:,loc(p,0,0)).*(p*(1+y)+F3);
    dby(:,loc(p,1,0)) = (p+F3y).*db(:,loc(p,0,0)) +...
        dby(:,loc(p,0,0)).*(p*(1+y)+F3);
    dbz(:,loc(p,1,0)) = F3z.*db(:,loc(p,0,0)) +...
        dbz(:,loc(p,0,0)).*(p*(1+y)+F3);
end

for p=0:k-2
    for q=1:k-p-1
        db(:,loc(p,q+1,0))=(a(2*p+1,0,q)*F4+b(2*p+1,0,q)*F5).*db(:,loc(p,q,0))...
            -c(2*p+1,0,q)*F5.^2.*db(:,loc(p,q-1,0));
        dbx(:,loc(p,q+1,0))=(a(2*p+1,0,q)*F4x+b(2*p+1,0,q)*F5x).*db(:,loc(p,q,0))...
            +(a(2*p+1,0,q)*F4+b(2*p+1,0,q)*F5).*dbx(:,loc(p,q,0))...
            -2*c(2*p+1,0,q)*F5.*F5x.*db(:,loc(p,q-1,0))...
            -c(2*p+1,0,q)*F5.^2.*dbx(:,loc(p,q-1,0));
        dby(:,loc(p,q+1,0))=(a(2*p+1,0,q)*F4y+b(2*p+1,0,q)*F5y).*db(:,loc(p,q,0))...
            +(a(2*p+1,0,q)*F4+b(2*p+1,0,q)*F5).*dby(:,loc(p,q,0))...
            -2*c(2*p+1,0,q)*F5.*F5y.*db(:,loc(p,q-1,0))...
            -c(2*p+1,0,q)*F5.^2.*dby(:,loc(p,q-1,0));
        dbz(:,loc(p,q+1,0))=(a(2*p+1,0,q)*F4z+b(2*p+1,0,q)*F5z).*db(:,loc(p,q,0))...
            +(a(2*p+1,0,q)*F4+b(2*p+1,0,q)*F5).*dbz(:,loc(p,q,0))...
            -2*c(2*p+1,0,q)*F5.*F5z.*db(:,loc(p,q-1,0))...
            -c(2*p+1,0,q)*F5.^2.*dbz(:,loc(p,q-1,0));
        
    end
end
    
for p=0:k-1
    for q=0:k-p-1
        db(:,loc(p,q,1)) =(1+p+q+(2+q+p)*z).*db(:,loc(p,q,0));
        dbx(:,loc(p,q,1)) =(1+p+q+(2+q+p)*z).*dbx(:,loc(p,q,0));
        dby(:,loc(p,q,1)) =(1+p+q+(2+q+p)*z).*dby(:,loc(p,q,0));
        dbz(:,loc(p,q,1)) =(1+p+q+(2+q+p)*z).*dbz(:,loc(p,q,0))...
            +(2+p+q)*db(:,loc(p,q,0));
    end
end

for p=0:k-2
    for q=0:k-p-2
        for r=1:k-p-q-1
            db(:,loc(p,q,r+1))=(a(2*p+2*q+2,0,r)*z+b(2*p+2*q+2,0,r)).*...
                db(:,loc(p,q,r))-c(2*p+2*q+2,0,r)*db(:,loc(p,q,r-1));
            dbx(:,loc(p,q,r+1))=(a(2*p+2*q+2,0,r)*z+b(2*p+2*q+2,0,r)).*...
                dbx(:,loc(p,q,r))-c(2*p+2*q+2,0,r)*dbx(:,loc(p,q,r-1));
            dby(:,loc(p,q,r+1))=(a(2*p+2*q+2,0,r)*z+b(2*p+2*q+2,0,r)).*...
                dby(:,loc(p,q,r))-c(2*p+2*q+2,0,r)*dby(:,loc(p,q,r-1));
            dbz(:,loc(p,q,r+1))=(a(2*p+2*q+2,0,r)*z+b(2*p+2*q+2,0,r)).*...
                dbz(:,loc(p,q,r))-c(2*p+2*q+2,0,r)*dbz(:,loc(p,q,r-1))...
                + a(2*p+2*q+2,0,r)*db(:,loc(p,q,r));
        end
    end
end



return

