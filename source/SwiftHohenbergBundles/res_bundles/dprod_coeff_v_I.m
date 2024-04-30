% assume indexing on v_coeff starts at -2 but v_{-1} = v_{-2} = 0.
function [coeff] = dprod_coeff_v_I(v_coeff, max_order, eigenvalues, growth_rate)
    lam_dot_alpha = lam_dot_alpha_mat(eigenvalues, max_order); 
    
    offset = 3;
    coeff = zeros(max_order + offset, max_order+offset); 
    for alpha = -2:max_order 
        for a1 = -2:alpha 
            a2 = alpha - a1; 
            if a2 > max_order
                continue 
            else 
                lam_a12 = lam_dot_alpha(a1 + offset, a2 + offset);
                
                const = (growth_rate + lam_a12)^3 - lam_a12 - growth_rate; 
                coeff(a1 + offset, a2 + offset) = const*v_coeff(a1 + offset, a2 + offset);
            end
            
        end
    end
end
