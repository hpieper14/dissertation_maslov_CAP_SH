offset = 3; 
w_coeff = zeros(3,3); 
w_coeff(-1 + offset, 0 + offset) = 1; 
w_coeff(0 + offset, -2 + offset) = 1; 
w_coeff(0+offset, -1 + offset) = 2; 
w_coeff(0 + offset, 0 + offset) = 3; 


params.eigenvalues.s = [-1, -1]; 
params.order = 0;

deriv_coeff = diff_w_coeff(params, w_coeff, 1); 


x_range = [0:.1:5]; 
[x_range, dw_vals] = diff_w(params, deriv_coeff, x_range);
[x_range, w_vals] = diff_w(params, w_coeff, x_range);

w = @(x) exp(2.*x) + 3.*exp(x) + 3; 
dw = @(x) 2.*exp(2.*x) + 3.*exp(x); 

true_w_vals = w(x_range); 
true_dw_vals = dw(x_range); 

% figure 
% tiledlayout(2,1)
% nexttile 
% hold on
% plot(x_range, dw_vals)
% plot(x_range, true_dw_vals)
% legend("with power series", "by_hand")
% title("$w'(x)$", Interpreter='latex')
% nexttile 
% hold on 
% plot(x_range, w_vals)
% plot(x_range, true_w_vals)
% legend("with power series", "by_hand")
% title("$w(x)$", Interpreter= "latex")

assert(sum(abs(w_vals - true_w_vals)) < 1e-10)
assert(sum(abs(dw_vals - true_dw_vals)) < 1e-10)

