% x is a 4 dim vector
function [vectors, values]= getJacEigs_toMerge(x, params)
    % compute the jacobian for SH
    Df0=JacSH_toMerge(x(1),params);
    [V1,D1]=eigs(Df0);
    % sort the eigenvalues and eigenvectors by increasing real part
    [d, ind]=sort(diag(real(D1)));
    D=D1(ind,ind);
    V=V1(:,ind);
    
    % scale the eigenvectors if necessary
    Vscale=params.scale*V;
    
    % save the eigenvalues and eigenvectors
    values.s=[D(1,1),D(2,2)];
    vectors.s=Vscale(:,1:2);
    values.u=[D(3,3),D(4,4)];
    vectors.u=Vscale(:,3:4);
    
    if (((real(values.s(1))<0) == 1) == 0) || (((real(values.s(2))<0) == 1) == 0)
        disp('There may be something wrong with the assignment of the asymptotic stable and unstable vectors')
    end
    
end
