% time should be in the form [min max]
% coeffs an mxn matrix with m-1 the order of the taylor expansion and n the
% number of components 
function [t,y] = taylorCoeffToFunction(coeffs, time)
    k = size(coeffs, 2);
    order = size(coeffs, 1)-1;
    
    t = time(1):.01:time(2);
    m = length(t);
    
    if isa(coeffs(1,1), 'intval') 
        y = intval(zeros(m,k));
    else
        y = zeros(m,k);
    end
    
    for i=1:k
        for j = 1:m
            coeffi = coeffs(:,i);
            y(j,i) = sum(coeffi.*(t(j).^(0:order))');
        end             
    end
end
