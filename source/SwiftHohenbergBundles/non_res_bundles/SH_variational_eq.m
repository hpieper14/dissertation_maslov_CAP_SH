function vals = SH_variational_eq(params, x_range, sol)
    N = max(size(x_range));
    vals = zeros(4, N); 
    SH_mat = SH_variational_mat(params, x_range);
    for i = 1:1:N
        this_sol = sol(i, :); 
        vals(:, i) = SH_mat(:,:,i)*this_sol';
    end
end
