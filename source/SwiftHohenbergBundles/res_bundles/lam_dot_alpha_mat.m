% inputs: lam = (lam1, lam2) 
% computes alpha \cdot \lam for |\alpha| \geq -2

function [coeff] = lam_dot_alpha_mat(lam, max_order)
    offset = 3;
    vec = -2:1:max_order; 
    rows = repmat(vec*lam(1)', [max_order + offset, 1])';
    cols = repmat((vec*lam(2)), [max_order + offset, 1]);

    coeff = rows + cols; 
  