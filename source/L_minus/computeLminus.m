function params = computeLminus(params, Lminus_cands) 
    r = (1+params.mu^2)^(1/2);
    upperbound = (7/(16*sqrt(r)))^(1/2)*(2*sqrt(2)*sqrt((1+2*sqrt(r))*(1/r^2 + 1)))^(-1);

    [~, values]= getBinfEigs(params);
    lam1 = values.u(1);
    lam2 = values.u(2); 
    
    shiftB = B_infinity(params) - real(lam1).*diag(ones(4,1));
    [V,~] = eigs(shiftB); 
    
    K = max(abs(V),[],'all')*max(abs(V^(-1)),[],'all');
    
    theta = (lam1 + lam2);
     
    % need to address the tail terms 
    norms = zeros(params.mfld.order, params.mfld.order);
    for i = 1:params.mfld.order 
        for j = 1:params.mfld.order
            norms(i,j) = vecnorm(params.unstable.coeffs(i,j,:));
        end
    end
    supQ = vecnorm(vecnorm(norms));
    
    C = supQ*(2*params.nu + 6*supQ);
    
    smallEnough = 0; 
    
    i = 1;
    while smallEnough == 0
        cand = Lminus_cands(i);
        tau = K*C/theta*exp(-theta*cand);
        
        if tau/(1-tau) < upperbound && tau < 1
            smallEnough = 1;
            disp('Found L_minus.')
            disp(Lminus_cands(i))
        else
            i = i + 1;
        end
    end
    params.Lminus = Lminus_cands(i); 
end
