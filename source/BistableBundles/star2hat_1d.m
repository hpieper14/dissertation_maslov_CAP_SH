% a, b are vectors of size (1, 2*order +1) with the coefficients arranged
% as (-order, order) 
function coeff = star2hat_1d(a,b, order)
    coeff = 0;
    for i = -order:order-2
        coeff = coeff + a(i + (order + 1))*b(order - i + (order + 1));
    end
end
