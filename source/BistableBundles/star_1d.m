% computes the terms of the a*b up to order m 
function coeffs = star_1d(a,b,m)
    assert(max(size(a)) >= m+1 & max(size(b)) >= m+1)
    coeffs = zeros(1, m + 1);
    for i = 0:m 
        coeff_i = 0;
        for j = 0:i
            coeff_i = coeff_i + a(j+1)*b(i-j+1);
        end
        coeffs(i+1) = coeff_i;
    end
end