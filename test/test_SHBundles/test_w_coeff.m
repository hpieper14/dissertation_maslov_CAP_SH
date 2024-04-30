x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/3;

order = 8; 

params.order = order; 

params.mfld.order = order; 

st_coeffs = getStBundleCoefficients(params);
v_coeff_pm = reshape(st_coeffs(1,:,:), [order + 1, order + 1]); 

[eigenvectors, eigenvalues] = getJacEigs_toMerge(0, params); 
st_mfld_coeffs = calc_proj_coeff(eigenvalues.s, eigenvectors.s, params);
mfld_deriv_coeff = st_mfld_deriv(params); 
v_coeff_dp = reshape(mfld_deriv_coeff(1,:,:), [order, order]); 

bottom = zeros(1, order); 
right = zeros(order + 1, 1); 
top = zeros(2, order + 1); 
left = zeros(order + 3, 2);
v_coeff_dp = [v_coeff_dp; bottom]; 
v_coeff_dp = [v_coeff_dp, right]; 
v_coeff_dp = [top; v_coeff_dp];
v_coeff_dp = [left, v_coeff_dp]; 



% pad v so it starts at index -2 
v_coeff_pm = pad_coeff_to_new_min_index(v_coeff_pm, 0, -2);

% drop last rows so dim of v matches that of w
v_coeff_pm = v_coeff_pm(1:end-3, 1:end-3); 
v_coeff_dp = v_coeff_dp(1:end-3, 1:end-3); 


w_coeff_pm = resonant_dprod_coeff(eigenvalues, v_coeff_pm, order - 3, - eigenvalues.s(1));
w_coeff_dp = resonant_dprod_coeff(eigenvalues, v_coeff_dp, order - 3, - eigenvalues.s(1));

x_range = [-15:.1:-12];
[x_range, vals_pm] = w_ODE(w_coeff_pm, v_coeff_pm, x_range, params);
[x_range, vals_dp] = w_ODE(w_coeff_dp, v_coeff_dp, x_range, params);