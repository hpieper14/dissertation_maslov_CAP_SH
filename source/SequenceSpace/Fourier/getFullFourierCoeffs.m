% returns a vector of coeffs (a_{-k},..., a_{-1}, a_0, a_1, .... a_k)
% order = n
function coeffs = getFullFourierCoeffs(sol,time,order,L)
    coeffs=zeros(1, 2*order+1);
    for n = -order:order
        F = sol.*exp(-1i*n*pi.*time./L);
        coeffs(n+order+1) = real(trapz(time, F)/(2*L));
    end       
  
    T=-L:.25:L;
    f = 0;
    for n = -order:order
        f = f+real(coeffs(n+order+1)*exp(1i*n*pi.*T./L));
    end
 
    
    figure
    plot(time, sol);
    hold on
    plot(T,f)
    legend('Sol via NF equation','Fourier approximation')

end