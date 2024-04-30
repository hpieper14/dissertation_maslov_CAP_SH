function val=testFourierDF(params,order, side)
    coeffs = rand(1,order+1);
    if side == 1
        approx_df = approx_fourier_DFfull(coeffs, order, params,side);
        DF  = DFFourierNBC(coeffs,order, params,side);
    else
        coeffs = [flip(coeffs), coeffs(2:end)];
        approx_df = approx_fourier_DFfull(coeffs, order, params,side);
        DF = DFFourierNDBC(coeffs,order, params); 
    end

    val=approx_df-DF;

    if max(val,[],'all')>1e-4
        disp('There might be something wrong with DF.')
    end

end