function df=approx_fourier_DFfull(coeffs, order, params,side)
    epsilon=1e-10;
    if side == 1
        df=zeros(order+1, order+1);
        for i = 1:order+1
            partial_deriv=approx_fourier_DFpartial(coeffs, i, order, params, epsilon,side);
            df(:,i)=partial_deriv';
        end
    else 
        df = zeros(2*order+1, 2*order+1);
        for i = -order:order 
            index = i + order+1; 
            partial_deriv=approx_fourier_DFpartial(coeffs, index, order, params, epsilon,side);
            df(:,index)=partial_deriv';
        end
    end
    
        
end
