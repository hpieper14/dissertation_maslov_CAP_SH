function [y1, y2] = ub_deriv_coords_x(params,ub_coeff, x_range)
    y1 = zeros(1, max(size(x_range))); 
    y2 = zeros(1, max(size(x_range)));
    u1 = ub_coeff(1,:); 
    u2 = ub_coeff(2,:);
    order = size(ub_coeff,2);



    sigma_0 = params.sigma_0; 
    for i = 1:max(size(x_range))
        x = x_range(i); 
        karg = [1:1:order].*(exp(x)*sigma_0).^[0:1:order-1];
        y1(i) = sigma_0*exp(x)*u1*karg';
        y2(i) = sigma_0*exp(x)*u2*karg';
    end
end
