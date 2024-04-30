function vals = bistable_variational_eq(params, x_range, sol)
    A = [0,1;1,0]; 
    B = @(t) [0,0; -6.*exp(-t)./(1+exp(-t)).^2, 0];
    
    N = max(size(x_range));
    vals = zeros(2, N); 

    for i = 1:N 
        uv = sol(:, i);
        vals(:,i) = B(x_range(i))*uv + A*uv;
    end
end
