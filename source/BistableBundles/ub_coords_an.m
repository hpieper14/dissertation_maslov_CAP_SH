function [y1, y2] = ub_coords_an(x_range)
    y1 = exp(x_range)./((1+exp(x_range)).^2); 

    y2 = -exp(x_range).*(-1+exp(x_range))./((1+exp(x_range)).^3);

end
