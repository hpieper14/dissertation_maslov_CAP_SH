function coords = sb_coords_roo(x_range, params,prod_coeff, ub_coeff)
    coords = zeros(2, max(size(x_range))); 
    order = params.bundle.order; 
    u_coeff = prod_coeff./[-2:1:order-1];
    u_coeff(3) = prod_coeff(3); 

    u_coeff = u_coeff(1:end-1);

    sigma_0 = params.sigma_0;
    for i = 1:max(size(x_range))
        x = x_range(i); 
        exp_vec = (exp(x)*sigma_0).^[1:1:order+1];
        unstab_sol = ub_coeff(1,:)*exp_vec'; 

        karg = [1:1:order+1].*(exp(x)*sigma_0).^[1:1:order+1];
        unstab_sol_deriv = ub_coeff(1,:)*karg';

        exp_vec = (exp(x)*sigma_0).^[-2:1:order-2];
        exp_vec(3) = x; 
        u_sol = u_coeff*exp_vec'; 

        u_sol_deriv = prod_coeff(1:end-1)*exp_vec';

        coords(1,i) = unstab_sol*u_sol;
        coords(2,i) = unstab_sol*u_sol_deriv + unstab_sol_deriv*u_sol;
    end
end 
