function sol = getFunctionFromFourierCoeffs(coeffs, L, order)
T=-L:.05:L;
f = 0;
for n = -order:order
    f = f+real(coeffs(n+order+1)*exp(1i*n*pi.*T./L));
end
sol(:,1) = T';
sol(:,2) = f';

end
