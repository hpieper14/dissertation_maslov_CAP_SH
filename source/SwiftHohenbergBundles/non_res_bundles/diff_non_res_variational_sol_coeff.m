function [coeffs] = diff_non_res_variational_sol_coeff(params, v_coeff, growth_rate, deriv_order)
    order = params.order;
    
    eigenvalues.s = params.eigenvalues.s; 

    offset = 1;
    coeffs = zeros(size(v_coeff));
    for alpha = 0:1:order
        for a1 = 0:1:alpha 
            a2 = alpha - a1;
            deriv_coeff = 0; 
            for j = 0:deriv_order
                deriv_coeff = deriv_coeff + nchoosek(deriv_order, j)*growth_rate^(deriv_order - j)...
                *dot(eigenvalues.s, [a1, a2])^j; 
            end
            coeffs(a1 + offset, a2 + offset) = v_coeff(a1 + offset, a2 + offset)*deriv_coeff;
        end
    end
end
