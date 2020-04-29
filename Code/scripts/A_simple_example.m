% Last modified: 4/12/2020
% A simple example on how to use HDG3D library to solve 
% Poisson's equation with Dirichlet boundary condition

% Before we do anything, let us add the HDG3D library
% Type "addpath('../')"

% We start by clearing any pre-existing variables and
% loading the domain and its discretizations
clear
load('../meshes/sixT3dDir');
% The above command loads a cubic domain - [0,1]^3, which is
% discretized by 6 tetrahedral meshes - {T1, T2, T3, T4, T5, T6}
% among which the mesh size decreases as the index increases

% We specify a mesh here
T = T3;

% We next define the polynomial degree, denoted by 'k'
% HDG3D library supports arbritray polynomial degrees
k = 2;

% Now we define the exact solution 'u'
% 'ux', 'uy', 'uz' are the three partial derivatives w.r.t. the axes
% 'uxx' is the x-derivative of 'ux', similar rules for 'uyy' and 'uzz'
u = @(x,y,z) sin(x.*y.*z);
ux = @(x,y,z) y.*z.*cos(x.*y.*z);
uy = @(x,y,z) x.*z.*cos(x.*y.*z);
uz = @(x,y,z) x.*y.*cos(x.*y.*z);
uxx = @(x,y,z) -y.^2.*z.^2.*sin(x.*y.*z);
uyy = @(x,y,z) -x.^2.*z.^2.*sin(x.*y.*z);
uzz = @(x,y,z) -x.^2.*y.^2.*sin(x.*y.*z);

% With 'uxx', 'uyy', and 'uzz', we can define the first input, which is
% the forcing term 'f = -uxx - uyy - uzz'
f = @(x,y,z) -uxx(x,y,z)-uyy(x,y,z)-uzz(x,y,z);

% Another input is the Dirichlet data, which we denote by 'uD'
uD = u;
% Note that here we simply use 'u' for 'uD', but actually any function 'f'
% is OK if the restriction of 'f' on the Dirichlet boundary is equal to 'u'

% Before we call the HDG solver, we need some extra variables
% 'formulas' stores the quadrature rules
% 'Nelts' stores the number of elements
% 'tau' stores the stabilization functions of the HDG method, type
% 'help createTau3d' to see more (don't forget to 'addpath' first!)
formulas=checkQuadrature3d(k,0);
Nelts=size(T.elements,1);
tau=createTau3d(Nelts,2);

% We are now ready to call the HDG3d solver
% For now, please ignore the use of 'onefun' and 'zfun'
% For more complicated PDEs (e.g. diffusion-reaction, Neumann B.C.)
% we will need to replace 'onefun' and 'zfun', but they are good for now
onefun = @(x,y,z) 1 + 0*x;
zfun = @(x,y,z) 0*x;
Uh=HDG3d(onefun,zfun,f,tau,T,k,formulas,uD,zfun,zfun,zfun);

% Finally, we are good to check the L^2 error
% Vary the mesh 'T' and the polynomial degree 'k' to see the changes
err = errorElem(T,u,Uh,k,formulas{1})

% Exercise: write a script to loop among the meshes and polynomial degrees
% to observe their effect on the errors. 
% What would be the convergence rates?
