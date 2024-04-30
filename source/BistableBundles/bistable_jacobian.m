function matrix = get_bistable_jacobian(b, x1, x2)
    matrix = [0,1; -3*b*(x1 + x1^2) + 1/2*b, 0];
end