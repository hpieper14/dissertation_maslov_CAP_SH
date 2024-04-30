function bundle_coeff = bundle_coeff_from_st_mfld_deriv(params, diff_var)
    assert((diff_var == 1) || (diff_var == 2))

    eigenvalues.s = params.eigenvalues.s; 
    eigenvectors.s = params.eigenvectors.s; 
    order = params.mfld.order; 

    mod_params = params; 
    mod_params.mfld.order = params.mfld.order + 1; 
    st_mfld_coeffs = calc_proj_coeff(eigenvalues.s, eigenvectors.s, mod_params);

    first_comp = reshape(st_mfld_coeffs(:, :, 1), [order + 2, order + 2]);
    vec = [1:1:order+1]; 
    mult_mat = ones(order+1, order+1); 
    mult_mat = mult_mat.*vec;

    if diff_var == 1
        mult_mat = mult_mat';
        first_comp = first_comp(2:end, 1:end-1);
    else 
        first_comp = first_comp(1:end-1, 2:end);
    end

    bundle_coeff = first_comp.*mult_mat;

    % lam_dot_alpha = lam_dot_alpha_mat(eigenvalues.s, params.order); 
    % lam_dot_alpha = lam_dot_alpha(3:end, 3:end);

    % st_mfld_coeffs(1,1,:) = zeros(1, 4); 

    % deriv_coeff1 = st_mfld_coeffs(:,:,1).*lam_dot_alpha';
    % deriv_coeff2 = st_mfld_coeffs(:,:,2).*lam_dot_alpha';
    % deriv_coeff3 = st_mfld_coeffs(:,:,3).*lam_dot_alpha';
    % deriv_coeff4 = st_mfld_coeffs(:,:,4).*lam_dot_alpha';

    % order = params.mfld.order; 
    
    % deriv_coeffs = zeros(4, order+1, order+1); 
    % deriv_coeffs(1, :, :) = reshape(deriv_coeff1, [order+1, order+1]); 
    % deriv_coeffs(2, :, :) = reshape(deriv_coeff2, [order+1, order+1]); 
    % deriv_coeffs(3, :, :) = reshape(deriv_coeff3, [order+1, order+1]); 
    % deriv_coeffs(4, :, :) = reshape(deriv_coeff4, [order+1, order+1]); 

    % deriv_coeffs(:, 1, 1) = zeros(4,1,1);
end
