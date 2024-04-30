function w_coeff = sb_product_deriv_coeff(ub_2dim, params)
    order = params.bundle.order;
    w_coeff = zeros(1, order + 2);
    w_coeff(1) = 1; 

    ub = ub_2dim(1,:);
    
    for k = 0:order
        sum_k = 0; 
        for k1 = -2:k-1 
            k2 = k - k1; 
            if k2 >=1 && k2 <= order
                sum_k = sum_k + (k1 + 2*k2)*w_coeff(k1+3)*ub(k2); 
            end
        end
        w_coeff(k+2) = -1/((k + 1)*ub(1))*sum_k;
    end

    w_coeff(6:end) = zeros(1, order - 3);
       
end
