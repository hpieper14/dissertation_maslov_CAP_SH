params.eigenvalues.s = [-1, -2]; 
params.order = 2; 

v_coeffs = [1, 2, 0; -1, 1, 0; 0, 0, 0]; 
dv_coeffs = diff_non_res_variational_sol_coeff(params, v_coeff, -1, 1);

x_range = 0:.1:5; 

[x_range, v_vals] = diff_non_res_variational_sol(params, v_coeffs, -1, x_range);
[x_range, dv_vals] = diff_non_res_variational_sol(params, dv_coeffs, -1, x_range);


v = @(x) exp(-x) + 2.*exp(-3.*x) - exp(-2.*x) + exp(-4.*x); 
dv = @(x) -exp(-x) - 6.*exp(-3.*x) +2.*exp(-2.*x) - 4.*exp(-4.*x); 

true_v_vals = v(x_range); 
true_dv_vals = dv(x_range);

figure 
tiledlayout(2,1)
nexttile 
hold on
plot(x_range, dv_vals)
plot(x_range, true_dv_vals)
legend("with power series", "by_hand")
title("$v'(x)$", Interpreter='latex')
nexttile 
hold on 
plot(x_range, v_vals)
plot(x_range, true_v_vals)
legend("with power series", "by_hand")
title("$v(x)$", Interpreter= "latex")

assert(sum(abs(v_vals - true_v_vals)) < 1e-10)
assert(sum(abs(dv_vals - true_dv_vals)) < 1e-10)
