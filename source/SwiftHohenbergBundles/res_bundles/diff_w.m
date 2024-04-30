function [x, vals] = diff_w(params, diff_w_coeff, x)
    order = params.order; 
    eigenvalues.s = params.eigenvalues.s;
    num_x = max(size(x)); 
    lam_dot_alpha = lam_dot_alpha_mat(eigenvalues.s, order);

    vals = zeros(1, num_x);
    for i = 1:1:num_x 
        v_sig_mat = exp(lam_dot_alpha.*x(i)).*diff_w_coeff;
        this_sum = sum(sum(v_sig_mat)); 
        vals(i) = this_sum;
    end
end
