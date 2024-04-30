function [coeff] = diff_w_coeff(params, w_coeff, deriv_order)
    deriv_coeff = (lam_dot_alpha_mat(params.eigenvalues.s, params.order)).^deriv_order;
    coeff = deriv_coeff.*w_coeff; 