function [y1, y2] = mfld_coords_an(x_range)
    y1 = 1./(1 + exp(-x_range));
    y2 = exp(-x_range)./(1 + exp(-x_range)).^2;
end
