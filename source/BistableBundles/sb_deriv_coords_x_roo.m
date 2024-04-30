function coords = sb_deriv_coords_x_roo(x_range, params, w_prime_coeff, ub_coeff)
    coords = zeros(2, max(size(x_range))); 
    order = params.bundle.order; 
    w_coeff = w_prime_coeff./[-2:1:order-1];
    w_coeff(3) = w_prime_coeff(3); 

    w_coeff = w_coeff(1:end-1);

    sigma_0 = params.sigma_0;
    for i = 1:max(size(x_range))
        x = x_range(i); 
        exp_vec = (exp(x)*sigma_0).^[1:1:order+1];
        unstab_sol = ub_coeff(1,:)*exp_vec'; 

        karg = [1:1:order+1].*(exp(x)*sigma_0).^[1:1:order+1];
        unstab_sol_deriv = ub_coeff(1,:)*karg';

        k2arg = ([1:1:order+1].^2).*(exp(x)*sigma_0).^[1:1:order+1];
        unstab_sol_2deriv = ub_coeff(1,:)*k2arg';

        exp_vec = (exp(x)*sigma_0).^[-2:1:order-2];
        exp_vec(3) = x; 
        w_sol = w_coeff*exp_vec'; 

        w_sol_deriv = w_prime_coeff(1:end-1)*exp_vec';
        w_sol_2deriv = w_prime_coeff(1:end-1)*(exp_vec.*[-2:1:order-2])';


        coords(1,i) = unstab_sol*w_sol_deriv + unstab_sol_deriv*w_sol;
        coords(2,i) = unstab_sol*w_sol_2deriv + 2*w_sol_deriv*unstab_sol_deriv ...
            + w_sol*unstab_sol_2deriv;

    end
end
