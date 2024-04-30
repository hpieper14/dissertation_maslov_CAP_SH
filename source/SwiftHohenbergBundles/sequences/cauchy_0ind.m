% a, b: structs with fields 
    % int order: max order 
    % array coeff: vector of coeffs
    % int index_lb: lowest order 
    % int dim: dimension of coeffs (i.e. dim of multi-index) 
% This paper looks promising for FFT implementation -- http://www1.cs.columbia.edu/~stratos/research/fft.pdf
function [a1a2] = cauchy_0ind(a1,a2,order)
    % assert dim(a1) = dim(a2)
    % assert order isn't too big -- cannot be the order of a1 or a2 because
    % we do not necessarily have enough terms to compute it 

    % can use FFT, but need to reference this post for indexing issues: 
    % https://www.mathworks.com/matlabcentral/answers/38066-difference-between-conv-ifft-fft-when-doing-convolution

    %fft_this = fftn(a1.coeff).*fftn(a2.coeff);
    

    all_coeff = convn(a1.coeff, a2.coeff); 
    if a1.dim == 1
        a1a2.coeff = all_coeff(1:order+1);
    elseif a1.dim == 2 
        trimmed_coeff = all_coeff(1:order+1, 1:order+1);
        a1a2.coeff = flip(tril(flip(trimmed_coeff)));
    else 
        disp("Multi-indices of dimension greater than 3 not currently supported.")

    end 
    
    a1a2.index_lb = a1.index_lb + a2.index_lb; 
    a1a2.order = order; 
    a1a2.dim = a1.dim; 
end 



