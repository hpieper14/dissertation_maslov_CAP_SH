
% Inputs: a,b - vectors of the same size 
% Outputs: (a*b)_k - scalar, where k+1 is the length of a,b and we are viewing a,b
% as the first k+1 components of infinite sequences 
function convolution = cauchy2(a,b,k) 
    if class(a) == 'intval'
        convolution=intval(0);
    else
        convolution = 0;
    end
    for j=0:k
        convolution=convolution+a(k-j+1)*b(j+1); 
    end
end