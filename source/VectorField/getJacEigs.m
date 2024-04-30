% x is a 4 dim vector
function [vectors, values]= getJacEigs(x, params)
    % compute the jacobian for SH
    Df0=JacSH(x(1),params.mu,params.nu);
    [V1,D1]=eigs(Df0);
    % sort the eigenvalues and eigenvectors by increasing real part
    [d, ind]=sort(diag(real(D1)));
    if real(D1(ind(1))) < real(D1(ind(end)))
        ind = flip(ind);
    end
    D=D1(ind,ind);
    V=V1(:,ind);
    
    % scale the eigenvectors if necessary
    Vscale=params.scale*V;
    
    % save the eigenvalues and eigenvectors
    values.u=[D(1,1),D(2,2)];
    vectors.u=Vscale(:,1:2);
    values.s=[D(3,3),D(4,4)];
    vectors.s=Vscale(:,3:4);
    
    if ((real(values.s(1))<0) == 1) == 0
        disp('There may be something wrong with the assignment of the asymptotic stable and unstable vectors')
    end
    
end
