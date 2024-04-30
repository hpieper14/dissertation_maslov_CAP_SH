function coords = ub_coords_pm(params)
    order = params.bundle.order;
    params = jacobian_eigenpairs(params);

    unstable_bundle_coeffs = ub_coeff(params);
    sigma_vec = 0:.05:1;
    p = size(sigma_vec, 2); 
    coords = zeros(p,2);

    for i = 1:p 
        sig = sigma_vec(i);
        p_i = zeros(2, 1);
        for n = 0:order 
            v_n = unstable_bundle_coeffs(:, n+1);
            p_i = p_i + v_n*sig^n; 
        end 
        coords(i,:) = sig*p_i;
    end
end
