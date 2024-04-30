function [vectors, values] = getBinfEigs(params)
    mat=B_infinity(params);
    [V1,D1]=eigs(mat);
    % sort the eigenvalues and eigenvectors by increasing real part
    [d, ind]=sort(real(diag(D1)));
    D=D1(ind,ind);
    V=V1(:,ind);
    
    % scale the eigenvectors if necessary

    
    % save the eigenvalues and eigenvectors
    values.s=[D(1,1),D(2,2)];
    vectors.s=V(:,1:2);
    values.u=[D(3,3),D(4,4)];
    vectors.u=V(:,3:4);
    
    if ((real(values.s(1))<0) == 1) == 0
        disp('There may be something wrong with the assignment of the asymptotic stable and unstable vectors')
    end
    
        
end