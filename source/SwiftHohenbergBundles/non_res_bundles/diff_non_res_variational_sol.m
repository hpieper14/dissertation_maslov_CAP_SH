function [x, vals] = diff_non_res_variational_sol(params, variational_deriv_coeff, growth_rate, x)
    order = params.order; 
    eigenvalues.s = params.eigenvalues.s;
    num_x = max(size(x)); 
    lam_dot_alpha = lam_dot_alpha_mat(eigenvalues.s, order);
    lam_dot_alpha = lam_dot_alpha(3:end, 3:end);

    vals = zeros(1, num_x);
    for i = 1:1:num_x 
        v_sig_mat = exp(lam_dot_alpha.*x(i)).*variational_deriv_coeff;
        this_sum = sum(sum(v_sig_mat)); 
        vals(i) = exp(growth_rate*x(i))*this_sum;
    end
end
