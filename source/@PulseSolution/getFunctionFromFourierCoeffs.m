% GETFUNCTIONFROMFOURIERCOEFFS  Computes function values from Fourier
% coefficients on either [-L,L] or [0,L]
%   sol = getFunctionFromFourierCoeffs(S,coeffs, interval_type)
%   sol = S.getFunctionFromFourierCoeffs(coeffs, interval_type)
% 
% returns a matrix of size 2 x N where the first column is the vector of
% times and the second column is the vector of function values. 
function sol = getFunctionFromFourierCoeffs(S,coeffs, interval_type)
order = S.fourier.order; 

if interval_type == "half"
    T=0:.05:S.time;
else
    T=-S.time:.05:S.time;
end

f = 0;
for n = -order:order
    f = f+real(coeffs(n+order+1)*exp(1i*n*pi.*T./S.time));
end

sol(:,1) = T';
sol(:,2) = f';

end
