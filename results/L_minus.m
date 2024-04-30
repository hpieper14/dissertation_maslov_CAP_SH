clear all; 
close all; 
% unstable 

disp("Calculating L minus for mu = 0.05 and nu = 1.6.")
params.nu = 1.6; 
params.mu = 0.05; 
params.lambda = 0; 


params.mfld.order = 15; 
[vectors, values] = getBinfEigs(params);
params.scale = 1/3; 

params.unstable.coeffs = calc_proj_coeff(values.u, params.scale*vectors.u,params);

Lminus_cands = [1:1:50];
params = computeLminus(params, Lminus_cands) ;


% stable 

disp("Calculating L minus for mu = 0.20 and nu = 1.6.")

params.nu = 1.6; 
params.mu = 0.2; 
params.lambda = 0; 


params.mfld.order = 15; 
[vectors, values] = getBinfEigs(params);
params.scale = 1/3; 

params.unstable.coeffs = calc_proj_coeff(values.u, params.scale*vectors.u,params);

Lminus_cands = [1:1:50];
params = computeLminus(params, Lminus_cands) ;