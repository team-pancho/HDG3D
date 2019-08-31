function tau=createTau3d(Nelts,option)

% tau=createTau3d(Nelts,option)
% 
% Input:
%    Nelts    : number of elements
%    option   : (1) tau = 1
%               (2) SFHDG on the first face
%               (3) SFHDG on a random face
% Output:
%    tau      : 4 x Nelts
%
% Last modified: March 14, 2013

switch option
    case 1
        tau=ones(4,Nelts);
    case 2
        tau=[ones(1,Nelts);...
             zeros(3,Nelts)];
    case 3
        where=ceil(4*rand(1,Nelts));
        tau=sparse(where,1:Nelts,1);
        tau=full(tau);
end
return


