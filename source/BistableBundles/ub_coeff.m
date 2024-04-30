function coeffs = ub_coeff(params)
    order = params.bundle.order;
    params = jacobian_eigenpairs(params); 
    Df0 = bistable_jacobian(params.b, 0,0);

    mflds.coeff.u = mfld_coeff(params, 'unstable'); 
    
    lam1 = params.evalues.unstable; 
    
    A = DFQbundle(params, mflds); 
    
    coeffs = zeros(2,order + 1); 
    coeffs(:,1) = [1;1]*params.bundle.scale;
    
    for i = 1:order 
        s_i = starhatMat(A, coeffs(:,1:i+1), i); 
    
        v_i = ((i*lam1+lam1)*eye(2) - Df0)^(-1)*s_i;
        coeffs(:,i+1) = v_i; 
        
    end
end
