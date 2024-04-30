function coeff = sb_coeff(params, prod_coeff, stab_coeff)
    % does this make sense? What does the order of these coefficients mean
    % when some of them have the form $\sigma^nx$ and some have the form $x
    % \sigma^nx$? 
    
    order = params.bundle.order; 
    
    % dim of coeff and order of computed stab bundle coeff do not coincide
    coeff = zeros(2, order + 1);

    for i = 1:order 
        v_i = 0;
        for k1 = -2:i-1 
            k2 = order - k1;
            if k2 >= 1 && k2 <= order
                v_i = v_i + prod_coeff(k1 + 3)*stab_coeff(k2 + 1);
            end
        end
    coeff(1, i) = v_i; 
    end

    % need to differentiate to obtain the second coordinate 



end
