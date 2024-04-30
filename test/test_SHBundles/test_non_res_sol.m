x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/3;

order = 8; 
params.order = order; 
params.mfld.order = order; 
[eigenvectors, eigenvalues] = getJacEigs_toMerge(0, params); 

params.eigenvalues.s = eigenvalues.s; 
params.eigenvectors.s = eigenvectors.s; 
params.eigenvalues.u = eigenvalues.u; 
params.eigenvectors.u = eigenvectors.u; 

st_coeffs = getStBundleCoefficients(params);

v_coeff_pm = reshape(st_coeffs(1,1,:,:), [order + 1, order + 1]); 
v_coeff_dp = bundle_coeff_from_st_mfld_deriv(params, 1);

growth_rate = eigenvalues.s(1);
x_range = .0001:.1:15; 

[x_range, u_dp] = scalar_sol_to_4D_sol(params, v_coeff_dp, growth_rate, x_range);
[x_range, u_pm] = scalar_sol_to_4D_sol(params, v_coeff_pm, growth_rate, x_range); 

u_dp = real(u_dp); 
u_pm = real(u_pm);

dv_coeff_dp = diff_non_res_variational_sol_coeff(params, v_coeff_dp, growth_rate, 1);
dv_coeff_pm = diff_non_res_variational_sol_coeff(params, v_coeff_pm, growth_rate, 1);

%% TODO Think about real vs complex valued quantities 
[new_x_range, dot_u_dp_c] = scalar_sol_to_4D_sol(params, dv_coeff_dp, growth_rate, x_range); 
[new_x_range, dot_u_pm] = scalar_sol_to_4D_sol(params, dv_coeff_dp, growth_rate, x_range); 

dot_u_dp = real(dot_u_dp_c); 
u_dp = real(u_dp);

params.mfld.st_coeffs = calc_proj_coeff(eigenvalues.s, eigenvectors.s, params); 
f_u_dp = SH_variational_eq(params, x_range, u_dp);
f_u_pm = SH_variational_eq(params, x_range, u_pm); 

should_be_0_dp = dot_u_dp - f_u_dp'; 
should_be_0_pm = dot_u_pm - f_u_pm'; 

disp("ERROR")
disp(vecnorm(should_be_0_pm))
disp(vecnorm(should_be_0_dp))


figure 
hold on 
plot(x_range, real(u_dp(:,1)))
plot(x_range, real(u_pm(:,1)))












