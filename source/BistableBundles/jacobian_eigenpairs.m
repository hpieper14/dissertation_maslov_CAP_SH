function params = jacobian_eigenpairs(params) 
    matrix = bistable_jacobian(params.b, 0, 0);
    [V1,D1] = eig(matrix); 
    [d, ind]=sort(diag(D1));
    D=D1(ind,ind);
    V=V1(:,ind);

    params.evectors.unstable = V(:,2);
    params.evectors.stable = V(:,1);
    params.evalues.unstable = D(2,2);
    params.evalues.stable = D(1,1);
end
