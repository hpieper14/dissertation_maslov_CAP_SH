% Inputs: a,b,c - vectors of the same size 
%         k - k entry in the cauchy product when using zero indexing
% Outputs: (a*b*c)_k - kth entry of the Cauchy product, is a scalar 
function convolution = cauchy3(a,b,c,k) 
    if class(a) == 'intval'
        convolution=intval(0);
    else
        convolution = 0;
    end
    
    for j=0:k
        for l=0:j
            convolution=convolution+a(k-j+1)*b(j-l+1)*c(l+1); 
        end
    end
end