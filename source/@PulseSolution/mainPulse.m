% MAINPULSE  Computes Fourier approximation of pulse solution and
% approximates the unstable eigenvalues associated to the pulse solution.
%   S = S.mainPulse()
%   S = mainPulse(S) 
% 
function S = mainPulse(S)
    % get pulse approximation from normal form. Determines where to
    % truncate spatial domain so $\varphi'(L) \approx 0$. 
    S = S.BKNormalForm4d_halfline();
    S = S.trimNFSol_halfline(); 
    S = BKNormalForm4d_halfline(S); 
    
    % refine Fourier coefficients with Newton's method
    S = S.Newton_halfline(); 

    % get Jacobian 
    DF = S.DFFourier(S.fourier.full_coeff_from_half_newton);
    
    % compute eigenvalues of Jacobian and save positive eigenvalues
    D = eig(DF);
    [i,j] = find(D>1e-10);
    vals = [];
    if length(i) > 0
        if min(size(D(i,j))) == 1
            vals = [vals, D(i,j)];
        else 
            A = D(i,j);
            vals = [vals, A(:, 1)];
        end
    end
    
    S.fourier.unstable_eigs = vals;
end