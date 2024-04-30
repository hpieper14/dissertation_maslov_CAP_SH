function y = chebcoeff_to_function(coeff)
    dom = -1:.05:1;
    n = max(size(coeff));
    m = max(size(dom));
    polys = zeros(m,n);
    
    polys(:,1) = ones(m,1);
    polys(:,2) = dom; 
    
    for i = 3:n
        Ti = 2.*dom'.*polys(:,i-1) - polys(:,i-2);
        polys(:,i) = Ti;
    end
    
    y = polys*coeff';
end
