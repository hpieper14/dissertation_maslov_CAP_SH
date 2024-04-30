% Non-Resonant Vector Bundles 
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

disp(['These are the first few coefficients for a stable vector bundle over' ...
    'over the stable manifold for the Swift-Hohemberg equation, calculated ' ...
    'by differentiating the parameterization of the stable manifold: '])

disp(v_coeff_dp(1:3, 1:3))

disp("These are the same coefficients calculated via the paramterization" + ...
    "method. We expect them to match the above coefficients.")
disp(v_coeff_pm(1:3, 1:3))
