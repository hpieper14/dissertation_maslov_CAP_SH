function prod = starMat(a,b,m,n)
    sizeA = size(a);
    sizeB = size(b); 
    if sizeA(end-2) ~= sizeB(end-3)
        msg = 'Error occurred. The inputs into starhatMat are not compatible for matrix multiplication.';
        error(msg)
    end
    
    prod = zeros(sizeA(1), sizeB(2));
    for i = 0:m
        for j = 0:n
            l = m-i;
            k = n-j;
            prod = prod+a(:,:,l+1, k+1)*b(:,:,i+1, j+1);
        end
    end
end
