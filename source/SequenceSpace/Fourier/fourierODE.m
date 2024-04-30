function fcn = fourierODE(a,order,params, side)
    if (size(a,1) == 1) == 0
        a = a';
    end
    
    length=max(size(a));
    if (order+1 == length) == 0
        a=[a,zeros(1,order-length)];
    end
      
    fcn=[];
    atripstar=tripconvvec(a,a,a,order,side);
    adubstar=dubconvvec(a,a,order,side);
    if side == 1
        for k=0:order
            coeff= -params.mu - (1-(k*pi/params.fourier.L)^2)^2; 
            fcn=[fcn,coeff*a(k+1) - atripstar(k+1) + params.nu*adubstar(k+1)];
        end  
    else
        for k = -order:order
            coeff= -params.mu - (1-(k*pi/params.fourier.L)^2)^2;
            index = k+order+1;
            fcn=[fcn,coeff*a(index) - atripstar(index) + params.nu*adubstar(index)];
        end 
    end
end