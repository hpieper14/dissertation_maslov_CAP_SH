function vals = bistable_variational_eq_pm(params, x_range, sol)
    A = [0,1;1,0]; 
    sigma_0 = params.sigma_0; 
    mfld_coeffs = params.mfld.coeffs; 
    
    N = max(size(x_range));
    vals = zeros(2, N); 

    for i = 1:N 
        x = x_range(i); 
        mfld_coord = mfld_coeffs(:,1)'*((exp(x)*sigma_0').^[0:1:params.mfld.order])';
        B = [0,0; 6*mfld_coord^2 - 6*mfld_coord, 0];
        uv = sol(:, i);
        vals(:,i) = B*uv + A*uv;
    end
end
