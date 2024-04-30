x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/3;

order = 6; 
params.order = order; 
params.mfld.order = order; 

[eigenvectors, eigenvalues] = getJacEigs_toMerge(0, params); 
scaled_evec = eigenvectors.s(:,1)*params.scale;
params.eigenvalues = eigenvalues; 
params.eigenvectors = eigenvectors; 

st_coeffs = getStBundleCoefficients(params);
v_coeff_dp = bundle_coeff_from_st_mfld_deriv(params, 1);
scaled_v_coeff_dp = v_coeff_dp; 


v_coeff_pm = reshape(st_coeffs(1,1,:,:), [order + 1, order + 1]); 


% plot norms 
flat_order_vec = 0:1:params.order; 
for i = 1:params.order
    flat_order_vec = [flat_order_vec, (flat_order_vec + i)];
end

%flat_order_vec = repmat(order_vec, 1, params.order+1);

coeff_norms_pm = abs(reshape(v_coeff_pm', 1, [])); 
index_to_keep_pm = find(coeff_norms_pm); 
x_pm = flat_order_vec(index_to_keep_pm); 
y_pm = coeff_norms_pm(index_to_keep_pm); 

coeff_norms_dp = abs(reshape(scaled_v_coeff_dp', 1, []));
index_to_keep_dp = find(coeff_norms_dp); 
x_dp = flat_order_vec(index_to_keep_dp); 
y_dp = coeff_norms_dp(index_to_keep_dp);


figure 
tiledlayout(2,1)
nexttile
plot(x_pm, log(y_pm), 'o')
xlabel('order')
ylabel('log norm')
title('Coefficient norms using IVB Invariance Equation')
nexttile 
plot(x_dp, log(y_dp), 'o')
xlabel('order')
ylabel('log norm')
title('Coefficient norms differentiating manifold power series')









% assert the two sets of coefficients are the same 
should_be_eigenvector = reshape(st_coeffs(:,1,1), [4 1]);

diff_norm = norm(v_coeff_pm - v_coeff_dp);
assert(isequal(scaled_evec, should_be_eigenvector))
assert(isequal(scaled_v_coeff_dp, v_coeff_pm))






