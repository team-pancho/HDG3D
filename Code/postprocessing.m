function uhstar=postprocessing(T,km,qhx,qhy,qhz,uh,k,formula)

% uhstar=postprocessing(T,km,qhx,qhy,qhz,uh,k,formula)
%
% Last modified: April 9, 2013


% Dimensions

d3plus=nchoosek(k+1+3,3);
d3    =nchoosek(k+3,3);
Nelts =size(T.elements,1);

% Matrices

S=Stiffness(T,k+1,formula);
S(1,:,:)=[];
S(:,1,:)=[];
[Cx,Cy,Cz]=VariableConv(T,km,km,km,k+1,formula);
Cx=permute(Cx,[2 1 3]);
Cy=permute(Cy,[2 1 3]);
Cz=permute(Cz,[2 1 3]);
Cx(:,d3+1:end,:)=[];
Cy(:,d3+1:end,:)=[];
Cz(:,d3+1:end,:)=[];
Cx(1,:,:)=[];
Cy(1,:,:)=[];
Cz(1,:,:)=[];

% Solution of local problems

uhstar=zeros(d3plus-1,Nelts);
parfor K=1:Nelts
    uhstar(:,K)=-S(:,:,K)\(Cx(:,:,K)*qhx(:,K)+...
                               Cy(:,:,K)*qhy(:,K)+...
                               Cz(:,:,K)*qhz(:,K));
end
uhstar=[uh(1,:); uhstar];

