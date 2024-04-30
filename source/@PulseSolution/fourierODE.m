% FOURIERODE  computes Fourier transform of the first order system
% evaluated at a, where a is the vector of Fourier coefficients
%   fcn = fourierODE(S, a) 
%   fcn = S.fourierODE(a)  
% 
% returns a vector
function fcn = fourierODE(S, a) 
    order = S.fourier.order; 
    params = S.vfParams;
    L = S.time;



    if (size(a,1) == 1) == 0
        a = a';
    end
    
    length=max(size(a));
    if (order+1 == length) == 0
        a=[a,zeros(1,order-length)];
    end
      
    fcn=[];
    side = 2;
    convolution = FourierSeries(order);
    atripstar=convolution.fourierConv3(a,a,a,order,side);
    adubstar=convolution.fourierConv2(a,a,order,side);

    for k = -order:order
        coeff= -params.mu - (1-(k*pi/L)^2)^2;
        index = k+order+1;
        fcn=[fcn,coeff*a(index) - atripstar(index) + params.nu*adubstar(index)];
    end 

end