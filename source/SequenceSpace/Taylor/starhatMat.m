% Written by Hannah Pieper. Last modified 11/15/21

% This file performs the Cauchy product for two doubly indexed truncated
% sequences a,b of the same order and returns the (m,n) element of the product. 
% The vectors a,b can have entries that are scalars or intervals.

% Inputs: 
%   a,b - doubly indexed truncated sequences. Formatted as
%   multi-dimensional arrays in which the last two dimensions give the
%   (i,j)th coefficient. These could be vector sequences or matrix
%   sequences. The dimensions of a,b excluding the last two give the
%   dimensions of the vector or matrix resp. 
%   m, n - nonnegative integers. 
% Outputs: 
%   prod - the (m,n) entry of (a^*b). 

% b has two columns of coefficients for order (l,k)
function prod = starhatMat(a,b,m,n)
    sizeA = size(a);

    prod = zeros(sizeA(1), 1); 
    for i = 0:m 
        for j = 0:n 
            l = m-i; 
            k = n-j; 
            prod = prod + a(:,:,i+1, j +1)*b(:, :, l + 1, k + 1);
        end
    end



    % prod = zeros(sizeA(1), 1);
    % for i = 0:m
    %     for j = 0:n-1
    %         l = m-i;
    %         k = n-j;
    %         prod = prod+a(:,:,l+1, k+1)*b(:,i+1, j+1);
    %     end
    % end
    % 
    % for i = 0:m-1
    %     j = n;
    %     l = m-i;
    %     k = n-j;
    %     prod = prod+a(:,:,l+1,k+1)*b(:,i+1, j+1);
    % end
end
