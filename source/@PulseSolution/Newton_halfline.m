% NEWTON_HALFLINE  Perform Newton's method on the pulse solution
% coefficients on domain [0,L]
%   S = S.Newton_halfline()
%   S = Newton)_halfline(S) 
% 
function S = Newton_halfline(S) 
    % Option to display output and use Jacobian
    options=optimset('Display','iter','Jacobian','on','MaxIter',10000);     
    
    % Sometimes need to tweak the scaling of the normal form solution for
    % convergence to a pulse 

    u = S.normalForm.sol(:,1); 
    [uout,fval] = fsolve(@(u) S.fourierODE_halfline(u),u,options);  
    coeffs = S.getFullFourierCoeffs(uout, S.normalForm.time);
    S.fourier.half_coeffs = coeffs; 
    
    % save full coefficients via reflection 
    full_uout = [flip(uout); uout];
    full_time = [-flip(S.normalForm.time); S.normalForm.time];
    S.fourier.full_coeff_from_half_newton = S.getFullFourierCoeffs(full_uout, full_time);
end
