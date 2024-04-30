function [t, y] = sb_product_coords_roo(params, t, ub_2dim)
    order = params.bundle.order;
    w_coeff = sb_product_deriv_coeff(ub_2dim, params);
    coords = zeros(2, max(size(t))); 
    order = max(size(w_coeff)); 
    int_vec = 1./[-2:1:(order-3)];
    u_coeff = w_coeff.*int_vec; 
    u_coeff(3) = w_coeff(3); 
    
    sigma_0 = params.sigma_0; 
    for i = 1:max(size(t))
        x = t(i);
        exp_vec = (sigma_0*exp(x)).^[-2:1:order-3];
        exp_vec(3) = x; 
        u_sol = u_coeff*exp_vec'; 
        coords(1, i) = u_sol;
    end

    y = coords(1,:);
end
