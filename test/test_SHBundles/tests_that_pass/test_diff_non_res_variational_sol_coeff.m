params.eigenvalues.s = [1, 2]; 
params.order = 2; 

v_coeff = [1, 2, 0; -1, 1, 0; 0, 0, 0]; 
dv_coeffs = diff_non_res_variational_sol_coeff(params, v_coeff, 1, 1);
dv_coeff_by_hand = [1, 6, 0; -2, 4, 0; 0, 0, 0]; 
assert(isequal(dv_coeffs, dv_coeff_by_hand))