function [v_inv_coeff] = v_inv(v_coeff, order)

    all_coeff = zeros(order + 1, order + 1); 
    v_00 = v_coeff(1,1);

    for alpha = 0:order
        for i = 0:alpha 
            j = alpha - i; 
            sub_vinv = all_coeff(1:i+1, 1:j+1); 
            sub_v = v_coeff(1:i+1, 1:j+1);
            all_coeff(i+1, j+1) = 1/v_00 * conv(sub_vinv, sub_v);
        end
    end
    
    trimmed_coeff = all_coeff(1:order+1, 1:order+1);
    v_inv_coeff = flip(tril(flip(trimmed_coeff)));

end
