% NEWTON  Performs Newton's method on the Fourier coefficients of the pulse
% solution
%   S = S.Newton()
%   S = Newton(S) 
% 
function S=Newton(S)
    s = S.fourier.full_coeff_from_half_newton;     
    k=0;
    while k < 100
        fcn=S.fourierODE(s);
        DF=S.DFFourier(s);
        if vecnorm(fcn)< S.fourier.tol 
            disp('Tolerence met for Newtons method.')
            break
        end
        s = s - (DF\(fcn'))';
        k=k+1;
    end
    S.fourier.full_coeff=s;
end