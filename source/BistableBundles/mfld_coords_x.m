function [y1, y2] = mfld_coords_x(params, x_range)
    y1 = zeros(1, max(size(x_range))); 
    y2 = zeros(1, max(size(x_range)));
    u1 = params.mfld.coeffs(:,1); 
    u2 = params.mfld.coeffs(:,2);
    order = params.mfld.order;

    sigma_0 = params.sigma_0;  
    for i = 1:max(size(x_range))
        x = x_range(i); 
        arg = (exp(x)*sigma_0).^[0:1:order];
        y1(i) = arg*u1; 
        y2(i) =arg*u2;
    end
end
