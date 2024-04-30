% FOURIERODE_HALFLINE  Computes Fourier approximation of pulse solution on
% the half line
% 
%   [F,J] = fourierODE_halfline(S, u)
%   [F,J] = S.fourierODE_halfline(u)
% 
function [F,J]=fourierODE_halfline(S, u)

    % construct differentiation matrices
    M = S.fourier.M;
    dt = 2*pi/M;
    M2 = M/2; 
    
    D2t = toeplitz([-pi^2/(3*dt^2)-1/6 ...
         .5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
    
    % rewrite matrix for 0..pi reflect
    semiD2t = zeros(M2+1);
    semiD2t(:,1) = D2t(M2:M,M2);
    semiD2t(:,2:M2) = D2t(M2:M,M2-1:-1:1)+D2t(M2:M,M2+1:M-1);
    semiD2t(:,M2+1)=D2t(M2:M,M);
    
    D2x  = (pi/S.time)^2*semiD2t; % Differentiation matrices on half line 
    
    N = M2+1;
    I = eye(N);
    LN = -D2x^2 - 2*D2x - I; 
    
    mu = S.vfParams.mu; 
    nu = S.vfParams.nu; 
    
    
   F = LN*u - mu*u + nu*u.^2 - u.^3;

    
   J = LN + diag(-mu + 2*nu*u - 3*u.^2);


