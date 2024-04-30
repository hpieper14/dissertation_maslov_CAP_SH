function dydt = nonautonomousODE(t,y,params, pulse)
    matrix = B_infinity(params);
    
    exponential = exp(-exponent.*t^2);
    
    exp_mat = [0,0,0,0;0,0,0,0;exponential,0,0,0; 0,0,0,0]; 
    
    matrix = matrix + exp_mat; 
    
    if (size(y,1) == 4) == 0
        y = y';
    end
    
    if (params.normalize == 1) == 1
        dydt = matrix*y - y.*dot(matrix*y, y);
    else 
        dydt = matrix*y;
    end
end