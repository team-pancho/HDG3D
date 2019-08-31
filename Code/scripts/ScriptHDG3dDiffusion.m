%% 
%Last modified: March 18, 2014
% Start by opening the parallel toolbox

% matlabpool open     % COMMENT IF NO PARALLEL TOOLBOX

% Close at the end

% Script Experiments with HDG3d for Reaction-Diffusion problems

clear 

k=input('Polynomial degree (<=3): ');
example=input('Input example (1->P1, 2->P2, 3->P3, 4->Cinfty): ');
ctmass=input('Coefficients: (1) variable; (0) constant: ');
domain=input('Domain: (1) chimney (2) cube (3) corner: ');
switch domain
    case 1
        load('../meshes/4Tchimney')     % Chimney-shaped domain - Dirichlet = horizontal
    case 2
        load('../meshes/sixT3dDir')     % (0,1)^3 all Dirichlet
    case 3
        load('../meshes/Corner')        % Fichera corner (uniform partition) - mixex BC
end
listT={T1,T2,T3,T4};

switch ctmass
    case 1
        kappa = @(x,y,z) 2+sin(x).*sin(y).*sin(z);
        kx = @(x,y,z) cos(x).*sin(y).*sin(z);
        ky = @(x,y,z) sin(x).*cos(y).*sin(z);
        kz = @(x,y,z) sin(x).*sin(y).*cos(z);
        c = @(x,y,z) 1.+0.5*(x.^2+y.^2+z.^2);
    case 0
        kappa = @(x,y,z) 2+0.*x;
        kx = @(x,y,z) 0.*x;
        ky = @(x,y,z) 0.*x;
        kz = @(x,y,z) 0.*x;
        c = @(x,y,z) 1+0.*x;
end

formulas=checkQuadrature3d(k,ctmass);

switch example
    case 1
        u = @(x,y,z) 2*x+3*y+z;
        ux = @(x,y,z) 2+0.*x;
        uy = @(x,y,z) 3+0.*x;
        uz = @(x,y,z) 1+0.*x;
        uxx = @(x,y,z) 0.*x;
        uyy = @(x,y,z) 0.*x;
        uzz = @(x,y,z) 0.*x;
    case 2
        u = @(x,y,z) x.^2+y.*x+z.*y+2*z.^2;
        ux = @(x,y,z) 2*x + y;
        uy = @(x,y,z) x + z;
        uz = @(x,y,z) y + 4*z;
        uxx = @(x,y,z) 2+0.*x;
        uyy = @(x,y,z) 0.*x;
        uzz = @(x,y,z) 4+0.*x;
    case 3
        u = @(x,y,z) 2*x.^2.*y+4*y.^2.*z+3*x.*z.^2;
        ux = @(x,y,z) 3*z.^2 + 4*x.*y;
        uy = @(x,y,z) 2*x.^2 + 8*y.*z;
        uz = @(x,y,z) 4*y.^2 + 6*x.*z;
        uxx = @(x,y,z) 4*y;
        uyy = @(x,y,z) 8*z;
        uzz = @(x,y,z) 6*x;
    case 4
        u = @(x,y,z) sin(x.*y.*z);
        ux = @(x,y,z) y.*z.*cos(x.*y.*z);
        uy = @(x,y,z) x.*z.*cos(x.*y.*z);
        uz = @(x,y,z) x.*y.*cos(x.*y.*z);
        uxx = @(x,y,z) -y.^2.*z.^2.*sin(x.*y.*z);
        uyy = @(x,y,z) -x.^2.*z.^2.*sin(x.*y.*z);
        uzz = @(x,y,z) -x.^2.*y.^2.*sin(x.*y.*z);
end

uD = u;
km = @(x,y,z) 1./kappa(x,y,z);
qx = @(x,y,z) -kappa(x,y,z).*ux(x,y,z);
qy = @(x,y,z) -kappa(x,y,z).*uy(x,y,z);
qz = @(x,y,z) -kappa(x,y,z).*uz(x,y,z);

f = @(x,y,z) -(kx(x,y,z).*ux(x,y,z)+kappa(x,y,z).*uxx(x,y,z))...
             -(ky(x,y,z).*uy(x,y,z)+kappa(x,y,z).*uyy(x,y,z))...
             -(kz(x,y,z).*uz(x,y,z)+kappa(x,y,z).*uzz(x,y,z))...
             +c(x,y,z).*u(x,y,z);

gx = @(x,y,z) -qx(x,y,z);         
gy = @(x,y,z) -qy(x,y,z);  
gz = @(x,y,z) -qz(x,y,z); 

a = @(x,y,z) 0.*x;

ErrorU=[];
ErrorQ=[];
ErrorUhat=[];
ErrorPu=[];
ErrorPuhat=[];
ErrorUstar=[];
h=[];

% Norms of unknowns for relative error

Tmax=listT{end};
Nelts=size(Tmax.elements,1);
d3=nchoosek(k+3,3);
normU=errorElem(Tmax,u,zeros(d3,Nelts),k,formulas{1});
normQ=errorElem(Tmax,qx,zeros(d3,Nelts),k,formulas{1})...
        +errorElem(Tmax,qy,zeros(d3,Nelts),k,formulas{1})...
        +errorElem(Tmax,qz,zeros(d3,Nelts),k,formulas{1});
    
for i=1:length(listT)
    T=listT{i};
    Nelts=size(T.elements,1);
    tau=createTau3d(Nelts,2);    
    [Uh,Qxh,Qyh,Qzh,Uhat]=HDG3d(km,c,f,tau,T,k,formulas,uD,gx,gy,gz);
    normUhat=errorFaces(T,u,zeros(size(Uhat)),k,formulas{4});
    
    % Errors with exact solution
    error_uhat=errorFaces(T,u,Uhat,k,formulas{4});
    error_q   =errorElem(T,qx,Qxh,k,formulas{1})...
               +errorElem(T,qy,Qyh,k,formulas{1})...
               +errorElem(T,qz,Qzh,k,formulas{1});
    error_u   =errorElem(T,u,Uh,k,formulas{1});
    ErrorQ=[ErrorQ error_q/normQ]; 
    ErrorU=[ErrorU error_u/normU];
    ErrorUhat=[ErrorUhat error_uhat/normUhat];
    
    % Errors with projections
    Puhat=L2projskeleton3d(u,T,k,formulas{3},formulas{4});
    [Pqx,Pqy,Pqz,Pu]=projectHDG3d(T,{qx,qy,qz,u},k,tau,formulas);
    error_Puhat=errorFaces(T,a,Uhat-Puhat,k,formulas{4});
    error_Pu=errorElem(T,a,Pu-Uh,k,formulas{1});
    ErrorPu=[ErrorPu error_Pu/normU];
    ErrorPuhat=[ErrorPuhat error_Puhat/normUhat];
    
    % Postprocessing
    Uhstar=postprocessing(T,km,Qxh,Qyh,Qzh,Uh,k,formulas{1});
    error_Ustar=errorElem(T,u,Uhstar,k+1,formulas{1});
    ErrorUstar=[ErrorUstar error_Ustar/normU];
       
    h=[h 1/(2^i)];
end

rateQ=log2(ErrorQ(1:end-1)./ErrorQ(2:end));
rateU=log2(ErrorU(1:end-1)./ErrorU(2:end));
rateUhat=log2(ErrorUhat(1:end-1)./ErrorUhat(2:end));
ratePu=log2(ErrorPu(1:end-1)./ErrorPu(2:end));
ratePuhat=log2(ErrorPuhat(1:end-1)./ErrorPuhat(2:end));
rateUstar=log2(ErrorUstar(1:end-1)./ErrorUstar(2:end));

format bank
disp('Rate |Q-Qh|   Rate |U-Uh|     Rate |U-Uhat|_h');
disp([rateQ' rateU' rateUhat']);
disp('Rate |\Pi u-Uh|  Rate |Pu-Uhat|_h  Rate |U-U*|')
disp([ratePu' ratePuhat' rateUstar']);

format shortE
disp('Error |Q-Qh|   Error |U-Uh|     Error |U-Uhat|_h');
disp([ErrorQ' ErrorU' ErrorUhat']);
disp('Error |\Pi u-Uh|  Error |Pu-Uhat|_h Error |U-U*|')
disp([ErrorPu' ErrorPuhat' ErrorUstar']);





%%

% matlabpool close        % COMMENT IF NO PARALLEL TOOLBOX
