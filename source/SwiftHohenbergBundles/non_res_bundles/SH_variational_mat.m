function vals = SH_variational_mat(params, x_range)
    A = [0,1,0,0; 0,0,1,0; 0,0,0,1; -1, 0, -2, 0]; 
    mfld_coeffs = params.mfld.st_coeffs; 

    p1 = reshape(mfld_coeffs(:,:,1), [params.order + 1, params.order + 1]); 

    N = max(size(x_range));
    vals = zeros(4, 4, N); 

    eigenvalues.s = params.eigenvalues.s;
    lam_dot_alpha = lam_dot_alpha_mat(eigenvalues.s, params.order); 
    lam_dot_alpha = lam_dot_alpha(3:end, 3:end);

    % COMPUTE MFLD COORD USING CODE!! THIS SHOULD BE REAL VALUED
    for i = 1:1:N
        sig_x_p1_mat = exp(lam_dot_alpha.*x_range(i)).*p1;
        mfld_coord = sum(sum(sig_x_p1_mat)); 
        non_lin = 2*params.nu*mfld_coord - 3*mfld_coord^2 - params.mu; 
        B = [0,0,0,0;0,0,0,0;0,0,0,0;non_lin, 0, 0, 0]; 
        vals(:, :, i) = (A + B);
    end