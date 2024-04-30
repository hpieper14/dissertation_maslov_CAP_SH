params.order = 2; 
params.eigenvalues.s = [-1, -2]; 

v_coeffs = [1, 2, 0; -1, 1, 0; 0, 0, 0]; 
v = @(x) exp(-x) + 2.*exp(-3.*x) - exp(-2.*x) + exp(-4.*x); 
dv = @(x) -exp(-x) - 6.*exp(-3.*x) + 2.*exp(-2.*x) - 4.*exp(-4.*x); 
d2v = @(x) exp(-x) + 18.*exp(-3.*x) - 4.*exp(-2.*x) + 16.*exp(-4.*x); 
d3v = @(x) -exp(-x) - 54.*exp(-3.*x) + 8.*exp(-2.*x) - 64.*exp(-4.*x); 


x_range = 0:.1:5;
[x_range, test_vals] = scalar_sol_to_4D_sol(params, v_coeffs, -1, x_range);

true_vals = [v(x_range)', dv(x_range)', d2v(x_range)', d3v(x_range)'];

is_0 = abs(true_vals - test_vals); 

assert(norm(is_0)< 1e-10)
