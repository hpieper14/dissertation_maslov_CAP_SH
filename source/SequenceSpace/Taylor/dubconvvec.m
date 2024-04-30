% Written by Hannah Pieper. Last modified 11/14/21

% This file performs the Cauchy product for two vectors a,b of the same
% order. The vectors a,b can have entries that are scalars or intervals.

% Inputs: 
%   side - 0, 1. Pass in 1 if the vectors a,b are one-sided meaning we have
%       assumed a_{-k} = a_{k} and we choose to work with a=(a_0,..., a_n). Pass in
%       0 if a=(a_{-n}, ... a_{-1}, a_0, a_1, ... a_n).
%   a,b - vectors. If a is one-sided, must be of size order+1. If a is two-sided, 
%       2*order+1. Can be column or row vectors.
%   order - order of the sequence

% Outputs: 
%   If the inputs are one-sided, the output is a row vector of size 2*order+1
%   containing the convolution terms of \tilde a = [a, zeros(1,order)]. If 
%   the inputs are two sided, then the output
%   is a row vector of size 4*order+1 containing the convolution terms of
%   \tilde a = [zeros(1,order), a, zeros(1,order)]. 
%
%   This is done because the lower order terms of the sequence a and b affect
%   higher order terms in the derivative of the convolution, so the extra 
%   terms may be useful in subsequent computations. 
function vec = dubconvvec(a,b,order, side)
    if size(a,2) == 1
        a = a';
    end
    
    if size(b,2) == 1
        b = b';
    end

    
    vec=[];
    if side == 1
       if (max(size(a)) == 2*order+1) == 0
        a = [a, zeros(1, order)];
        b = [b, zeros(1, order)];
       end
        
        for k=0:order
            dubsum=0;
            for i=-order:order
                l=k-i;
                dubsum=dubsum+a(abs(i)+1)*b(abs(l)+1);
            end
            vec=[vec,dubsum];
        end
    else 
        if (max(size(a)) == 4*order+1) == 0
            a = [zeros(1, order), a, zeros(1, order)];
            b = [zeros(1, order), b, zeros(1, order)];
        end
        for k=-order:order
            dubsum=0;
            for i=-order:order
                l=k-i;
                
                m = i+2*order+1;
                n = l+2*order+1;
                
                dubsum=dubsum+a(m)*b(n);
            end
            vec=[vec,dubsum];
        end
    end
    
end