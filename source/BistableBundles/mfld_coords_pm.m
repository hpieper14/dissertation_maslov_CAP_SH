function coords = mfld_coords_pm(params)
    order = params.bundle.order;
    params = jacobian_eigenpairs(params); 
    mflds.coeff.u = mfld_coeff(params, 'unstable'); 

    sigma_vec = 0:.05:1;
    p = size(sigma_vec, 2);
    coords = zeros(p, 2);
    
    for i = 1:p 
        sig = sigma_vec(i);
        p_i = zeros(2, 1); 
        for n = 0:order
            q_n = mflds.coeff.u(n+1, :);
            p_i = p_i + q_n'*sig^n; 
        end
        coords(i, :) = p_i;
    end
end
