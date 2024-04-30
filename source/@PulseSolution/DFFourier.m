% DFFOURIER  computes the Jacobian for the Fourier transform of the first
% order system evaluated at a, where a is the vector of Fourier
% coefficients
%   mat = DFFourier(S, a)
%   mat = S.DFFourier(a) 
% 
% returns a matrix
function mat = DFFourier(S, a)
    order = S.fourier.order;
    params = S.vfParams; 
    
    % when computing DF_k/Da_j, we differentiate with respect to a_j and a_{-j}
    mat=zeros(2*order+1, 2*order+1);
    
    if (size(a,1) == 1) == 0
        a = a';
    end
    
    side = 2;
    
    a=[zeros(1,order), a, zeros(1,order)];
    convolution = FourierSeries(order);
    duba=convolution.fourierConv2(a,a, 2*order, side);
    
    L = S.time;
    mu = params.mu;
    nu = params.nu;
    
    % lth row corresponds to F_l
    % kth column corresponds to d/da_k
    
    
    for l = -order:order 
        matrowLind = l  + order + 1; 
        for k = -order:order
            matcolKind = k + order + 1; 
            aLKind = l - k + 2*order+1;
            mat(matrowLind, matcolKind) = -3*(duba(aLKind)) ...
                    + 2*nu*(a(aLKind));
            if k == l
                mat(matrowLind,matcolKind) = mat(matrowLind,matcolKind)-mu-(1-(k*pi/L)^2)^2;
            end
        end
end
        
