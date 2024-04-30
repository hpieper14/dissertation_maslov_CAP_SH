% params -- params.b, params.evectors.unstable/stable,
% params.evalues.unstable/stable, params.fp = 0
function coeff = mfld_coeff(params, stability)
    if strcmp(stability,'unstable')
        eigenvectors = params.evectors.unstable;
        eigenvalues = params.evalues.unstable;
    else
        eigenvectors = params.evectors.stable;
        eigenvalues = params.evalues.stable;
    end
    
    order=params.mfld.order;
    coeff=zeros(order+1,2);
    e1=params.mfld.scale*eigenvectors(:,1);
    lam1=eigenvalues(1);
    
    Df0=bistable_jacobian(params.b,0,0); 
    % the zeroth order coefficient is the equilibrium, corresponding to 0.
    coeff(1,:) = params.fp;
    coeff(2,:) = [1,1]*params.mfld.scale;
   for j=2:order 
       A_j = (Df0 - j*lam1*eye(2))^(-1);
       % can include the j+1 entry because this is the starhat product
       a = coeff(1:j+1, 1);
       aa_j = starhat_1d(a,a,j);
       aaa_j = dubstarhat_1d(a,a,a,j); 
       B = [0, params.b*3/2*aa_j - params.b*aaa_j];
       p_j = A_j*B' ;
       coeff(j+1, :) = p_j;
   end
end