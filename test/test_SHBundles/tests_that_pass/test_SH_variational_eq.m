% We linearize about the function f(x) = 1+2exp(-2x)-exp(-x)+exp(-3x)
params.order = 2; 
params.eigenvalues.s  = [-1, -2]; 
params.mu = .05; 
params.nu = 1.6;

mfld_coeffs = [1, 2, 0; -1, 1, 0; 0, 0, 0]; 

params.mfld.st_coeffs = zeros(4, params.order + 1, params.order +1); 
params.mfld.st_coeffs(1, :, :) = mfld_coeffs;

mfld = @(x) 1 + 2.*exp(-2.*x) - exp(-x) + exp(-3.*x); 

x_range = 0:.1:5; 
test_vals = SH_variational_mat(params, x_range);

num = max(size(x_range));
true_vals = zeros(4,4,num);
norms = zeros(1, num);
for i = 1:1:num 
    v_x = mfld(x_range(i)); 
    nonlinearity = -1 + params.nu*v_x^2 - v_x^3 - params.mu*v_x;
    true = [0,1,0,0; 
                0,0,1,0;
                0,0,0,1;
                nonlinearity,  0, -2, 0];
    true_vals(:,:,i) = true; 
    norms(i) = sum(sum(abs(true - test_vals(:,:,i))));
end 

assert(sum(norms) < 1e-10)



