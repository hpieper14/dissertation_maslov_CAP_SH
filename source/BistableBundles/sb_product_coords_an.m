function [t,y] = sb_product_coords_an(t)
    y = (1/2.*exp(-2.*t).*(-1-8.*exp(t)+8.*exp(3.*t) + exp(4.*t))+6.*t);
end
