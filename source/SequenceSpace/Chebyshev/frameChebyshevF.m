% the first four columns of cheb coeffs correspond to the pulse, the second
% four columns correspond to the basis solutions
function F = frameChebyshevF(coeffs, params, ic, L)
    order = size(coeffs,1)-1;
    
    a1 = coeffs(:,1); 
    a2 = coeffs(:,2); 
    a3 = coeffs(:,3); 
    a4 = coeffs(:,4); 
    
    a5 = coeffs(:,5);
    a6 = coeffs(:,6); 
    a7 = coeffs(:,7);
    a8 = coeffs(:,8); 
    
 
    % chebstar functions use differnt order convention
    a1a1 = chebstar2(a1,a1,order+1)';
    a1a1a1 = chebstar3(a1,a1,a1,order+1)';
    a1a1a5 = chebstar3(a1,a1,a5,order+1)'; 
    a1a5 = chebstar2(a1,a5,order+1)';
   
    
    c1 = [a2;0];
    c2 = [a3;0];
    c3 = [a4;0]; 
    c4 = [-2*a3-(params.mu + 1)*a1 + params.nu*a1a1 - a1a1a1;0];
    
    c5 = [a8;0];
    c6 = [a7-2*a8;0];
    test = (params.lambda - 1 - params.mu)*a5 +2*params.nu*a1a5 - 3*a1a1a5;
    c7 = [test; 0];
    c8 = [a6;0];
    
    c = [c1,c2,c3,c4,c5,c6,c7,c8];
    
    F = zeros(order+1, 8); 
    
    ones = (-1).^(1:order);
    
    for i = 1:8
        F(1,i) = coeffs(1,i) + 2*sum(coeffs(2:end,i).*ones') - ic(i);
    end
    
    
    for k = 1:order 
        for i = 1:8
        F(k+1, i) = 2*coeffs(k,i) - L/k*(c(k,i)- c(k+2,i));
        end
    end
   
end
