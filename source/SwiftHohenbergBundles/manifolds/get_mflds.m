% params must contain mu, nu, scale, mflds.order
function params = get_mflds(params)
    % get scaled eigenvectors 
    [vectors, values]= getJacEigs(0, params);
    
    params.stable.coeffs=calc_proj_coeff(values.s, vectors.s, params);
    params.unstable.coeffs=calc_proj_coeff(values.u, vectors.u, params);
    
    
end
