% A is matrix valued, size (2,2,m) and b is vector valued size (2,m) 
% computes the term (A*b)_m only 
function prod = starhatMat(A,b,m)
    sizeA = size(A);
    sizeB = size(b); 
    if sizeA(end-1) ~= sizeB(end-1)
        msg = 'Error occurred. The inputs into starhatMat are not compatible for matrix multiplication.';
        error(msg)
    end
    
%     prod = zeros(sizeA(1), m+1);
%     for i = 0:m
%         coeff_i = zeros(2,1);
%         for j = 0:i
%             coeff_i = coeff_i + A(:,:,j+1)*b(:,i-j+1);
%         end
%         prod(:, i+1) = coeff_i;
%     end

    prod = zeros(2, 1); 
    for j = 1:m 
        prod = prod + A(:,:,j+1)*b(:, m-j+1); 
    end
end