% sol_vec should be a Nx5 dimensional vector with the first column being
% time and the 2-5 columns the spatial coordinates.
function cheb_coeff = get_cheb_coeffs(sol_vec, params)
    M=params.cheb.order;
    step=1:1:M;
    
    % must call chebfun on a column vector
    x1=chebfun(sol_vec(:,1),'trunc', M);
    x2=chebfun(sol_vec(:,2),'trunc', M);
    x3=chebfun(sol_vec(:,3),'trunc', M);
    x4=chebfun(sol_vec(:,4),'trunc', M);
    
    a1=chebcoeffs(x1);
    a2=chebcoeffs(x2);
    a3=chebcoeffs(x3);
    a4=chebcoeffs(x4);
    
    cheb_coeff=[a1, a2, a3, a4];
    
    
    figure
    tiledlayout(4,1)
    nexttile
    plot(x1)
    nexttile
    plot(x2)
    nexttile
    plot(x3)
    nexttile
    plot(x4)
    title('Cheb Solution (Pre-Newton Method)')
    
    figure
    tiledlayout(4,1)
    nexttile
    plot(step,a1, 'o')
    nexttile
    plot(step,a2, 'o')
    nexttile
    plot(step,a3,'o')
    nexttile
    plot(step,a4,'o')
    title('Chebyshev Coefficients')
    
    
end
