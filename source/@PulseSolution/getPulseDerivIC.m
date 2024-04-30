% GETPULSEDERIVIC  computes the initial conditions to the solution of the
% first order system constructed from the derivative of the pulse at time
% t_0. 
%   S = S.getPulseDerivIC(t_0)
%   S = getPulseDerivIC(S, t_0) 
% 
function S = getPulseDerivIC(S, t_0)
    L = S.time; 
    if t_0 > L 
        error('Must set initial time to be in [-L, L]')
    end 

    order = S.fourier.order;
    
    chainRuleTerms = -1i*pi*[-order:order]/L;
    ak = S.fourier.full_coeff_from_half_newton; 
    expTerms = ak.*exp(-1i*pi*[-order:order]*t_0/L);
    

    derivCoeffs = zeros(4, 1); 
    for i = 1:4 
        derivCoeffs(i, 1) = (chainRuleTerms.^i) * expTerms'; 
    end

    symplecticCoords = [derivCoeffs(1), derivCoeffs(3), derivCoeffs(4) + 2*derivCoeffs(2), derivCoeffs(2)];

    S.fourier.pulseDeriv.coords = real(symplecticCoords); 
    S.fourier.pulseDeriv.t_0 = t_0; 
end 
