function [x, vals] = w_ODE(w_coeff, v_coeff, x, params)
    % TODO: should change the diff_w functions to compute the coefficients inside of the function
    dw3_coeff = diff_w_coeff(params, w_coeff, 3); 
    dw3 = diff_w(params, dw3_coeff,x); 
    dw2_coeff = diff_w_coeff(params, w_coeff, 2); 
    dw2 = diff_w(params, dw2_coeff, x); 
    dw1_coeff = diff_w_coeff(params, w_coeff, 1); 
    dw1 = diff_w(params, dw1_coeff, x); 
    w = diff_w(params, w_coeff, x); 

    [~, values] = getJacEigs_toMerge(0, params);
    growth_rate = values.s(1);
    dv3_coeff = diff_non_res_variational_sol_coeff(params, v_coeff, growth_rate, 3);
    dv3 = diff_non_res_variational_sol(params, dv3_coeff, growth_rate, x); 
    dv2_coeff = diff_non_res_variational_sol_coeff(params, v_coeff, growth_rate, 2); 
    dv2 = diff_non_res_variational_sol(params, dv2_coeff, growth_rate, x); 
    dv1_coeff = diff_non_res_variational_sol_coeff(params, v_coeff, growth_rate, 1); 
    dv1 = diff_non_res_variational_sol(params, dv1_coeff, growth_rate, x); 
    v = diff_non_res_variational_sol(params, v_coeff, growth_rate, x);

    vals = -dw3.*v - 4*dw2.*dv1 + dw1.*(6.*dv2 - 2.*v) - 4*w.*(dv3 + dv1);
end
