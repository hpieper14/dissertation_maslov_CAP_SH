% Written by Hannah Pieper. Last modified 11/15/21

% This file performs the Cauchy product for two doubly indexed truncated
% sequences a,b of the same order and returns the (m,n) element of the product. 
% The vectors a,b can have entries that are scalars or intervals.

% Inputs: 
%   a,b - doubly indexed truncated sequences. Formatted as matrices
%       so that a_{i,j} is in a(i+1, j+1). 
%   m, n - nonnegative integers. 
% Outputs: 
%   prod - the (m,n) entry of (a^*b). 
function prod = starVec(a,b,m,n)
    testIndex = size(a);
    if isa(a(testIndex), 'intval')
        sum = intval(0);
    else
        sum=0;
    end
    
    for i=0:m
        for j=0:n
            sum=sum+a(m-i+1,n-j+1)*b(i+1,j+1);
        end
    end
    prod=sum;
end