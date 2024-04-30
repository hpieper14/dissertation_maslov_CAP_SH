offset = 3; 
w_coeff = zeros(3,3); 
w_coeff(-1 + offset, 0 + offset) = 1; 
w_coeff(0 + offset, -2 + offset) = 1; 
w_coeff(0+offset, -1 + offset) = 2; 
w_coeff(0 + offset, 0 + offset) = 3; 


params.eigenvalues.s = [-1, -1]; 
params.order = 0;

deriv_coeff = diff_w_coeff(params, w_coeff, 1); 
by_hand_deriv_coeff = [0,0,0; 0,0,1; 2, 2,0]; 

assert(isequal(deriv_coeff, by_hand_deriv_coeff))