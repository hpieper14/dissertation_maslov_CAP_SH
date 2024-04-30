function val=approx_fourier_DFpartial(coeffs, var, order, params, epsilon, side)
    perturb=zeros(1,max(size(coeffs)));
    perturb(var)=epsilon;
    perturb_coeffs=coeffs+perturb;
    
    val=(fourierODE(perturb_coeffs,order,params,side)-fourierODE(coeffs,order,params,side))/epsilon;
end
