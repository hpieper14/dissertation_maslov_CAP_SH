function matrix = B_infinity(params)
    matrix = [0,0,0,1;0,0,1,-2;params.lambda-1-params.mu,0,0,0;0,1,0,0];
end