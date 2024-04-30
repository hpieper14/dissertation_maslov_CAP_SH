% FOURIERCONV3  Computes the Cauchy product for three vectors a, b, c. 
%   vec = fourierConv2(a, b, order, side)
%
% Inputs: 
%   side - 0, 1. Pass in 1 if the vectors a,b are one-sided meaning we have
%       assumed a_{-k} = a_{k} and we choose to work with a=(a_0,..., a_n). Pass in
%       0 if a=(a_{-n}, ... a_{-1}, a_0, a_1, ... a_n).
%   a,b,c - vectors. If a is one-sided, must be of size order+1. If a is two-sided, 
%       2*order+1. Can be column or row vectors.
%   order - order of the sequence

% Outputs: 
%   If the inputs are one-sided, the output is a row vector of size 2*order+1
%   containing the convolution terms of \tilde a = [a, zeros(1,order)]. If 
%   the inputs are two sided, then the output
%   is a row vector of size 4*order+1 containing the convolution terms of
%   \tilde a = [zeros(1,order), a, zeros(1,order)]. 
function [vec] = fourierConv3(a,b,c,order, side)
    vec = [];
    if side == 1
        if (max(size(a)) == 3*order+1) == 0
            a = [a, zeros(1, 2*order)];
            b = [b, zeros(1, 2*order)];
            c = [c, zeros(1, 2*order)];
        end
        for k = 0:order 
            sum = 0;
            for k1 = -order:order
                for k2 = -order:order
                    k3 = k-k1-k2;
                    sum = sum+a(abs(k1)+1)*b(abs(k2)+1)*c(abs(k3)+1);
                end
            end
            vec = [vec, sum];
        end
    else
         if (max(size(a)) == 6*order+1) == 0
            a = [zeros(1,2*order), a, zeros(1, 2*order)];
            b = [zeros(1,2*order), b, zeros(1, 2*order)];
            c = [zeros(1,2*order), c, zeros(1, 2*order)];
         end
         for k = -order:order
             sum = 0;
             for k1 = -order:order
                 for k2 = -order:order
                     k3 = k - k1 - k2;
                     m = k1+3*order+1;
                     n = k2+3*order+1;
                     p = k3+3*order+1;
                     sum = sum+a(m)*b(n)*c(p);
                 end
            end
            vec = [vec, sum];
         end
    end
                     
        
    
end