function [vecs,eigs]=asym_un_vecs(params)
    [V1, D1] = eig(B_infinity(params));
    [d, ind]=sort(real(diag(D1)));
    D=D1(ind,ind);
    V=V1(:,ind);
    if (real(D(3)) < 0) || (real(D(4))< 0)
        disp('There might be a problem with the eigenvalues and vectors')
        disp(D)
    end
    eigs=[D(3,3),D(4,4)];
    vecs=V(:,3:4);
end