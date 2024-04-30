% GET4DIMIC  Computes the initial condition of the solution to the 4
% dimensional first order system corresponding to the pulse solution from
% its Fourier coefficients and saves the initial conditions to the
% ConjugatePoints object C
%
%   C = get4DimIC(C, S)
%   C = C.get4DimIC(S)
function C = get4DimIC(C, S)
    old_L = S.time;
    new_L = C.conjPts.L;

   % fourierCoeffs = S.fourier.full_coeff; 
    fourierCoeffs = S.fourier.full_coeff_from_half_newton; 
    order = S.fourier.order;
    derivative_coeffs = zeros(3, max(size(fourierCoeffs)));

    for k = -order:order
        index = k + order + 1; 
        for j = 1:3
            derivative_coeffs(j, index) = (1i*k*pi/old_L)^j*fourierCoeffs(index); 
        end
    end
  
    sol = S.getFunctionFromFourierCoeffs(fourierCoeffs, "full");
    solution_vec = sol(:,2); 

    for j = 1:3
        sol = S.getFunctionFromFourierCoeffs(derivative_coeffs(j,:), "full"); 
        solution_vec = [solution_vec, sol(:,2)]; 
    end

    time_vec = sol(:,1); 
    index = find(diff(sign(time_vec + new_L))); 
    if max(size(index)) > 1
        index = index(1);
    end

    C.Euminus.pulse_ic = solution_vec(index, :);
end