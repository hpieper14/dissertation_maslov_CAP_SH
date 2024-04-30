function dydt = asymptotic_ODE(t,y,params)
    matrix = B_infinity(params);
    
    if (size(y,1) == 4) == 0
        y = y';
    end
    
    if (params.normalize == 1) == 1
        dydt = matrix*y - y.*dot(matrix*y, y);
    else 
        dydt = matrix*y;
    end
end
