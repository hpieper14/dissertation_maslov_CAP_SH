% BINF  Calculates asymptotic coefficient matrix for given vector field
% parameter values. 
% 
%   matrix = Binf(C)
%   matrix = C.Binf()
function matrix = Binf(C)
    matrix = [0,0,0,1;...
        0,0,1,-2;...
        C.vfParams.lambda-1-C.vfParams.mu,0,0,0;...
        0,1,0,0];
end