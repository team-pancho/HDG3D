function formulas=checkQuadrature3d(k,constantmass)

%formulas=checkQuadrature3d(k,constantmass)
%Imput: 
%              k: degree of polynomials
%   constantmass: 1 constant coefficients in mass matrices
%                 0 non-costant
%Output:
%       formulas{1} : 3d quadrature formula for errors and var coeff
%       formulas{2} : 3d quadrature formula for constant coefficients
%       formulas{3} : 2d quadrature formula
%       formulas{4} : 2d quadrature formula for errors
%Last update: March 14, 2013

% degrees = {3k,2k,2k,2k+2} if constant mass = 0
%           {2k,2k,2k,2k+2} if constant mass = 1
% for k=0, take {2,0,0,2}

TablesQuadForm3d
TablesQuadForm
switch constantmass
    case 1
        switch k
            case 0
                formulas={tetra5,tetra1,matrix0,matrix4};
            case 1
                formulas={tetra3,tetra3,matrix4,matrix4};
            case 2
                formulas={tetra7,tetra7,matrix9,matrix9};
            case 3
                formulas={tetra9,tetra9,matrix11,matrix11};
        end
    case 0
        switch k
            case 0 
                formulas={tetra5,tetra1,matrix0,matrix4};
            case 1
                formulas={tetra5,tetra3,matrix4,matrix4};
            case 2
                formulas={tetra7,tetra7,matrix9,matrix9};
            case 3
                formulas={tetra9,tetra9,matrix11,matrix11};
        end
end
formulas{3}(:,4)=formulas{3}(:,4)/2;
formulas{4}(:,4)=formulas{4}(:,4)/2;
return