function [x_range, vals] = scalar_sol_to_4D_sol(params, coeff, growth_rate, x_range)
    disp(params)
    disp(growth_rate)
    % get derivative coefficients
    d_coeff = diff_non_res_variational_sol_coeff(params, coeff, growth_rate, 1); 
    d2_coeff = diff_non_res_variational_sol_coeff(params, coeff, growth_rate, 2);
    d3_coeff = diff_non_res_variational_sol_coeff(params, coeff, growth_rate, 3); 
    
    % compute functions and derivatives on the domain x_range using the
    % coefficients
    [x,sol] = diff_non_res_variational_sol(params, coeff, growth_rate, x_range); 
    [x,d_sol] = diff_non_res_variational_sol(params, d_coeff, growth_rate, x_range); 
    [x,d2_sol] = diff_non_res_variational_sol(params, d2_coeff, growth_rate, x_range); 
    [x,d3_sol] = diff_non_res_variational_sol(params, d3_coeff, growth_rate, x_range); 
    vals = [sol', d_sol', d2_sol', d3_sol']; 
end
