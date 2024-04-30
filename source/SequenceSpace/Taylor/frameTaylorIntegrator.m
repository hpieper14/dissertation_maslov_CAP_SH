function coeffs =  frameTaylorIntegrator(ic, order, params)
    if isa(ic(5), 'intval')
        coeffs = intval(zeros(order+1,8));
    else 
        coeffs = zeros(order+1, 8);
    end
    coeffs(1,:) = ic; 
        
    for i = 1:order 
        % compute coefficients for the pulse
        a1 = coeffs(1:i,1);
        coeffs(i+1,1) = coeffs(i,2);
        coeffs(i+1,2) = coeffs(i,3);
        coeffs(i+1,3) = coeffs(i,4);
        coeffs(i+1,4) = -2*coeffs(i,3) - (params.mu+1)*a1(i) ...
            + params.nu*cauchy2(a1, a1, i-1) ...
            - cauchy3(a1,a1,a1,i-1);
        coeffs(i+1,5) = coeffs(i,8);
        coeffs(i+1,6) = coeffs(i,7) - 2*coeffs(i,8);
        coeffs(i+1,7) = (-1-params.mu)*coeffs(i,5) + 2*params.nu*cauchy2(a1,coeffs(1:i,5),i-1) ...
                        -3*cauchy3(a1,a1,coeffs(1:i,5),i-1);
        coeffs(i+1,8) = coeffs(i,6);
        
        coeffs(i+1,:) = coeffs(i+1,:)/i;
    end   
end
