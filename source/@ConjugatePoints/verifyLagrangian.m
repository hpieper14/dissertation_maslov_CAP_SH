% VERIFYLAGRANGIAN Verifies that a given frame matrix is Lagrangian.
%
%   isLagrangian = verifyLagrangian(C, frame) 
%   isLagrangian = C.verifyLagrangian(frame) 
function isLagrangian = verifyLagrangian(C, frame) 
    if (size(frame, 1) == 2*size(frame, 2)) == 0
        disp('This frame may not be of size 2n x n') 
    end
    
    n = size(frame, 2); 
    
    X = frame(1:n, 1:n);
    Y = frame(n+1:end, 1:n);
    
    val = X'*Y - Y'*X;
    
    if (max(val) < 1e-13)
        isLagrangian = 1;
        disp('This subspace is Lagrangian.')
    else
        isLagrangian = 0;
        disp('This subspace is not Lagrangian.')
    end

end
