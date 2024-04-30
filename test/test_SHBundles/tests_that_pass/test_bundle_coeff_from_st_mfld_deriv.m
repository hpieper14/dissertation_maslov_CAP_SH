params.eigenvalues.stable = [-1, -2]; 
mfld_first_comp_coeffs = [1,-2,3,0;
                        2,3,0,0;
                        4,0,0,0;
                        0,0,0,0]; 

true_bundle_coeffs_1 = [2,3,0; 8,0,0; 0,0,0];

order = 2; 
vec = [1:1:order+1]; 
mult_mat = ones(order+1, order+1); 
mult_mat = mult_mat.*vec;
mult_mat = mult_mat';

test_bundle_coeff = mfld_first_comp_coeffs(2:end, 1:end-1).*mult_mat;

assert(isequal(test_bundle_coeff, true_bundle_coeffs_1))